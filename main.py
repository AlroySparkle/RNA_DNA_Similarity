import math
import csv
import os
import plotly.graph_objects as go
import tkinter as tk
from tkinter import filedialog
import tabulate as tabulate
import random
import matplotlib.pyplot as plt
import gc
import statistics
from datetime import datetime


def get_entropy(seq, window_size):
    entropy = []
    tsallis_entropy = []
    renyi_entropy = []
    window_open = 0
    window_close = window_open + int(window_size) - 1
    while window_close < len(seq):
        a_count = 0
        g_count = 0
        c_count = 0
        t_count = 0
        window_seq = seq[window_open:window_close]
        for char in window_seq:
            if char == 'A':
                a_count += 1
            else:
                if char == 'G':
                    g_count += 1
                else:
                    if char == 'C':
                        c_count += 1
                    else:
                        if char == 'T':
                            t_count += 1
        entropy.append(round(shannon(a_count, c_count, g_count, t_count, window_size), 2))
        renyi_entropy.append(round(renyi(a_count, c_count, g_count, t_count, 0.5), 2))
        tsallis_entropy.append(round(tsallis(a_count, c_count, g_count, t_count, 0.5), 2))
        window_open += 1
        window_close += 1
    return entropy, renyi_entropy, tsallis_entropy


def shannon(a_count, c_count, g_count, t_count, window_size):
    return -(a_count / window_size * math.log((a_count if a_count != 0 else 1) / window_size, 2)) - \
           (g_count / window_size * math.log((g_count if g_count != 0 else 1) / window_size, 2)) - \
           (t_count / window_size * math.log((t_count if t_count != 0 else 1) / window_size, 2)) - \
           (c_count / window_size * math.log((c_count if c_count != 0 else 1) / window_size, 2))


def renyi(a_count, c_count, g_count, t_count, double):
    window_size = a_count + c_count + g_count + t_count
    entropy = pow(a_count / window_size, double)
    entropy = entropy + pow(c_count / window_size, double)
    entropy = entropy + pow(g_count / window_size, double)
    entropy = entropy + pow(t_count / window_size, double)
    entropy = math.log2(entropy)
    return (1 / (1 - double)) * entropy


def tsallis(a_count, c_count, g_count, t_count, double):
    window_size = a_count + c_count + g_count + t_count
    entropy = pow(a_count / window_size, double)
    entropy = entropy + pow(c_count / window_size, double)
    entropy = entropy + pow(g_count / window_size, double)
    entropy = entropy + pow(t_count / window_size, double)
    entropy = 1 - entropy
    return (1 / (double - 1)) * entropy


def open_file():
    return str(filedialog.askopenfilenames()).replace('(', '').replace(')', '') \
        .replace("'", "").strip().split(',')


def file_save(file_name, shannons, renyis, tsallisses, save_directory):
    os.chdir(save_directory)
    with open(file_name + '.csv', 'w', newline='') as csv_file:
        fieldnames = ['id', 'shannon', 'renyi', 'tsallis']
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        for index in range(0, len(shannons)):
            writer.writerow({'id': index, 'shannon': shannons[index],
                             "renyi": renyis[index], "tsallis": tsallisses[index]})


def graph(entropies, renyis, tsalliss, e, r, t, names, seq_end):
    fig = go.Figure()
    length = []
    seq_start = 0
    for i in range(seq_start, seq_end):
        length.append(i)
    for index in range(len(entropies)):
        if e:
            fig.add_trace(go.Scatter(x=length, y=entropies[index], mode='lines', name=names[index] + ' shannon'))
        if r:
            fig.add_trace(go.Scatter(x=length, y=renyis[index], mode='lines', name=names[index] + ' renyi'))
        if t:
            fig.add_trace(go.Scatter(x=length, y=tsalliss[index], mode='lines', name=names[index] + ' tsallis'))
    fig.update_layout(title='result')
    fig.show()


def goParallel(list1, list2, window):
    if len(list1) < window or len(list2) < window:
        return
    temp = 0
    for i in range(window):
        temp += abs(list1[i] - list2[i])
    temp2 = temp
    last = min(len(list1), len(list2))
    for j in range(window, last):
        temp -= abs(list1[j - window] - list2[j - window])
        temp += abs(list1[j] - list2[j])
        temp2 = min(temp2, temp)
    return temp2


def walk(list1, list2, window):
    result = -1
    for incline in range(0, len(list1) + len(list2) - 2 * window + 1):
        if incline < len(list1) - window + 1:
            temp = goParallel(list1[len(list1) - window - incline:len(list1)], list2, window)
        else:
            temp = goParallel(list1, list2[incline - len(list1) + window:len(list2)], window)
        if not (temp is None or result > temp or result < -1):
            result = temp
    return result


def similarity(entropies, seq_length, names, file_name):
    outer_list = []
    inner_list = []
    for region1 in range(len(entropies)):
        for region2 in range(len(entropies)):
            if region1 == region2:
                result = 0.0
            elif region2 < region1:
                result = outer_list[region2][region1 + 1]
            else:
                result = walk(entropies[region1], entropies[region2], seq_length) / seq_length
            if len(inner_list) == 0:
                name = names.pop(0)
                inner_list.append(name)
                names.append(name)
            inner_list.append(round(result, 2))
        outer_list.append(inner_list)
        inner_list = []
    saving_files_list = [["names"] + names]
    for index in outer_list:
        saving_files_list.append(index)
    similarity_save(f=file_name, data_list=saving_files_list)
    # TODO make it return entire list later for specific test
    return outer_list


def select(choice):
    choice = choice.lower()
    return choice[0] == 'y', choice[1] == 'y', choice[2] == 'y'


def similarity_save(f, data_list):
    with open(f + ".csv", 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(data_list)
    f.close()


def make_bin(list1, bins):
    list1 = sorted(list1)
    max_val = max(list1)
    min_val = min(list1)
    bin_list = [0] * bins
    recent_index = 0
    bin_size = (max_val-min_val)/bins
    recent_size = bin_size+min_val
    index = 0
    while index < len(list1):
        if list1[index] > recent_size:
            recent_index += 1
            recent_size += bin_size
            continue
        bin_list[recent_index] += 1
        index += 1
    return bin_list


def main():
    gc.enable()
    window_size = int(input("window size: "))
    print("select entropy by write y or n in xyz")
    print("x shannon")
    print("y renyi")
    print("x tsallis")
    shannon_list = []
    renyi_list = []
    tsallis_list = []
    names = []
    shortest = -1
    while True:
        selected_entropies = input("selected: ")
        if len(selected_entropies) == 3:
            break
        print("input should be 3 letters only")
    e, r, t = select(selected_entropies)
    root = tk.Tk()
    root.withdraw()
    openFile = open_file()
    save = filedialog.askdirectory(title="save directory")
    print("open files: ")
    for f in openFile:
        filename = open(f.strip(), 'r')
        name = f.split('/')[-1].strip(".fa")
        if not f.endswith(".fa") and not f.endswith(".fasta"):
            continue
        seq = ''
        for x in filename:
            if x[0] == '>':
                name += "_" + x.split('|')[1]
                continue
            else:
                seq += x.rstrip()
        filename.close()
        if len(seq) < 50 or len(seq) - window_size <= window_size:
            continue
        if len(seq) < shortest or shortest < 0:
            shortest = len(seq)
        entropies, renyis, tsalliss = get_entropy(seq, window_size)
        shannon_list.append(entropies)
        renyi_list.append(renyis)
        tsallis_list.append(tsalliss)
        names.append(name.strip(".fa"))
    #   file_save(name, entropies, renyis, tsalliss, save)
    os.chdir(save)
    # graph(shannon_list, renyi_list, tsallis_list, e, r, t, names, selected_start, shortest-window_size)
    print("similarity: ")
    if e:
        test = similarity(shannon_list, 72, names, "shannon_similarity")
        counter = 1
        print("bootstrapping+bin")
        for index in test:
            mean_value = []
            data = index[1:-1]
            for i in range(1000):
                list3 = []
                for i in range(0, len(data)):
                    list3.append(random.choice(data))
                mean_value.append(sum(list3) / len(list3))
            mean = statistics.mean(mean_value)
            standard_deviation = statistics.stdev(mean_value)
            print("bootstrap", index, round(mean-standard_deviation*2,3), round(mean+standard_deviation*2,3))
            mean = statistics.mean(index[1:-1])
            print("-----------------")
        plt.show()
    # TODO return the buttom lines to get results for tsallis and renyi
    # if t:
    #    print("tsallis table:")
    #    similarity(tsallis_list, window_size, names, "tsallis_similarity")
    # if r:
    #    print("renyi table:")
    #    similarity(renyi_list, window_size, names, "renyi_similarity")


if __name__ == "__main__":
    main()
