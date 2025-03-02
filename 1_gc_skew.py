#!/usr/bin/env python3

"""
script for calculating gc skew

Chris Brown
ctb@berkeley.edu

modificado por Jerson Alexander García-Zea
jagarcia@eafit.edu.co

# python gc_skew.py -f *
# python gc_skew.py -f -l -w -s --single --no-plot

# posicionarse dentro del directorio con los fasta
# genera un output llamad 'total_gcskew.txt' con los ori - ter

# please cite how: Brown, C. T., Olm, M. R., Thomas, B. C., & Bandfiel, J. F. (2016). Measurement of bacterial
# replication rates in microbial communities. Nature biotechnology, 34(12), 1256-1263. http://doi.org/10.1038/nbt.3074
# "format APA"

"""

# python modules
import os
import sys
import argparse
import numpy as np
from scipy import signal
from itertools import cycle, product

# plotting modules
from matplotlib import use as mplUse
mplUse('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['pdf.fonttype'] = 42
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

# ctb
from ctbBio.fasta import iterate_fasta as parse_fasta

# Función para plotear los resultados
def plot_two(title, subtitle, A, B, labels, legend, vert = False):
    fig, ax1 = plt.subplots()
    colors = ['0.75', 'b', 'r', 'c', 'y', 'm', 'k', 'g']
    a_colors = cycle(colors)
    b_colors = cycle(colors[::-1])
    a_label = cycle(legend[0])
    b_label = cycle(legend[1])

    for a in A:
        x, y = a
        ax1.set_ylabel(labels[0], labelpad = 3)
        ax1.set_xlabel(labels[-1])
        ax1.plot(x, y, c = next(a_colors), marker = 'o', ms = 4, label = next(a_label))

    if vert is not False:
        for i in vert:
            x, c = i
            ax1.axvline(x = x, c = c, label = next(a_label), linewidth = 2)

    ax2 = ax1.twinx()
    for b in B:
        x, y = b
        ax2.set_ylabel(labels[1], labelpad = 8)
        ax2.plot(x, y, c = next(b_colors), linewidth = 2, label = next(b_label))

    xmin = min([min(i[0]) for i in A] + [min(i[0]) for i in B])
    xmax = max([max(i[0]) for i in A] + [max(i[0]) for i in B])
    ax2.set_xlim(xmin, xmax)

    plt.suptitle(title, fontsize = 16)
    plt.title(subtitle, fontsize = 10)

    ax1.legend(loc = 'upper left', bbox_to_anchor=(0.55, -0.125), prop = {'size':8}, framealpha = 0.0)
    plt.legend(loc = 'upper right', bbox_to_anchor=(0.45, -0.125), prop = {'size':8}, framealpha = 0.0)

    pdf = PdfPages('%s.pdf' % title.replace(' ', '_'))
    pdf.savefig(bbox_inches = 'tight')
    plt.close()
    pdf.close()

def check_peaks(peaks, length):
    closest, farthest = int(length * float(0.45)), int(length * float(0.55))
    pairs = []
    for pair in list(product(*peaks)):
        tr, pk = sorted(list(pair), key = lambda x: x[1], reverse = False)
        a = (tr[0] - pk[0]) % length
        b = (pk[0] - tr[0]) % length
        pt = abs(tr[1] - pk[1])
        if (a <= farthest and a >= closest) or (b <=farthest and b >= closest):
            pairs.append([pt, tr, pk])
    if len(pairs) == 0:
        return [False, False]
    pt, tr, pk = sorted(pairs, reverse = True)[0]
    return [tr[0], pk[0]]

def find_ori_ter(c_skew, length):
    c_skew_min = signal.argrelextrema(np.asarray(c_skew[1]), np.less, order = 1)[0].tolist()
    c_skew_max = signal.argrelextrema(np.asarray(c_skew[1]), np.greater, order = 1)[0].tolist()

    if len(c_skew_min) == 0 or len(c_skew_min) == 0:
        return [False, False]
    else:
        c_skew_min = [[c_skew[0][i], c_skew[1][i]] for i in c_skew_min]
        c_skew_max = [[c_skew[0][i], c_skew[1][i]] for i in c_skew_max]
        ori, ter = check_peaks([c_skew_min, c_skew_max], length)
    return ori, ter

def gc_skew(name, length, seq, window, slide, plot_skew):
    replacements = {'G':1, 'C':-1, 'A':0, 'T':0, 'N':0}
    gmc = []
    for base in seq:
        try:
            gmc.append(replacements[base])
        except:
            gmc.append(0)

    gpc = [abs(i) for i in gmc]
    weights = np.ones(window)/window
    gmc = [[i, c] for i, c in enumerate(signal.fftconvolve(gmc, weights, 'same').tolist())]
    gpc = [[i, c] for i, c in enumerate(signal.fftconvolve(gpc, weights, 'same').tolist())]

    skew = [[], []]
    c_skew = [[], []]
    cs = 0
    for i, m in gmc[0::slide]:
        p = gpc[i][1]
        if p == 0:
            gcs = 0
        else:
            gcs = m/p
        cs += gcs
        skew[0].append(i)
        c_skew[0].append(i)
        skew[1].append(gcs)
        c_skew[1].append(cs)
    ori, ter = find_ori_ter(c_skew, length)

    if plot_skew is True:
        title = '%s GC Skew' % (name)
        subtitle = '(window = %s, slide = %s)' % (window, slide)
        labels = ['GC Skew', 'Cumulative GC Skew', 'Position on Genome (bp)']
        N = int(len(skew[0])/1000)
        if N != 0:
            skew = [skew[0][0::N], skew[1][0::N]]
        if ori is False:
            plot_two(title, subtitle, [skew], [c_skew], labels, [[labels[0]], [labels[1]]])
        else:
            plot_two(title, subtitle, [skew], [c_skew], labels, [[labels[0], 'Ori:%s' % ('{:,}'.format(ori)), 'Ter:%s' % ('{:,}'.format(ter))], [labels[1]]], vert = [(ori, 'r'), (ter, 'b')])
    return ori, ter, skew, c_skew

def parse_genomes(fastas, single):
    """
    Generador para parsear archivos fasta
    Si single es True, combina las secuencias en un archivo multifasta
    """
    if single is True:
        for genome in fastas:
            sequence = []
            try:
                for seq in parse_fasta(genome):
                    if len(seq) < 2:
                        print(f"Error al procesar {genome.name}: secuencia inesperada o formato incorrecto.", file=sys.stderr)
                        continue
                    sequence.extend(list(seq[1].upper()))
                if not sequence:
                    print(f"No se encontró ninguna secuencia válida en el archivo {genome.name}.", file=sys.stderr)
                    continue
                yield (genome.name.rsplit('.', 1)[0], len(sequence), sequence)
            except Exception as e:
                print(f"Error al leer el archivo {genome.name}: {e}", file=sys.stderr)
                continue
    else:
        for genome in fastas:
            try:
                for seq in parse_fasta(genome):
                    if len(seq) < 2:
                        print(f"Error al procesar {genome.name}: secuencia inesperada o formato incorrecto.", file=sys.stderr)
                        continue
                    ID = seq[0].split('>', 1)[1].split()[0]
                    yield (ID, len(seq[1]), list(seq[1].upper()))
            except Exception as e:
                print(f"Error al leer el archivo {genome.name}: {e}", file=sys.stderr)
                continue


def open_files(files):
    if files is None:
        return files
    if files[0] == '-':
        return (sys.stdin)
    return (open(i) for i in files)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = '# calculate gc skew and find Ori and Ter of replication')
    parser.add_argument('-f', nargs = '*', action = 'store', required = True, help = 'fasta(s)')
    parser.add_argument('-l', default = False, type = int, help = 'minimum contig length (default = 10 x window)')
    parser.add_argument('-w', default = 1000, type = int, help = 'window length (default = 1000)')
    parser.add_argument('-s', default = 10, type = int, help = 'slide length (default = 10)')
    parser.add_argument('--single', action = 'store_true', help = 'combine multi-fasta sequences into single genome')
    parser.add_argument('--no-plot', action = 'store_false', help = 'do not generate plots, print GC Skew to stdout')
    args = vars(parser.parse_args())

    fastas = open_files(args['f'])
    single, plot_skew = args['single'], args['no_plot']
    window, slide = args['w'], args['s']
    min_len = args['l']
    if min_len is False:
        min_len = 10 * window

    # Abrir el archivo de salida
    with open('total_gcskew.txt', 'w') as outfile:
        for name, length, seq in parse_genomes(fastas, single):
            if length < min_len:
                continue
            ori, ter, skew, c_skew = gc_skew(name, length, seq, window, slide, plot_skew)

            # Escribir en el archivo de salida
            outfile.write(f"{name};{ori};{ter}\n")

            if plot_skew is False:
                continue

