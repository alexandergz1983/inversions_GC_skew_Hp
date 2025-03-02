# Required inputs:
# 1. Mochiview file with added columns, AA sequence and orientation
#        >>>   <genomename>_HOCD.mochi
# 2. Fasta file.
#        >>>   <genomename>.fasta
# 3. Ter coordinates file
#        >>>  "ter_coords.txt"
#
# >>> Genome names must match:
#     - fasta file name (the name in the file can be different)
#    - mochiview file
#
# >>> Program will find and analyze EVERY mochiview/fasta file in the working directory.

import os
from time import gmtime, strftime

#################################################################################################
def filenames():
    '''Collects the names of all mochiview files'''
    all_files = os.listdir()
    all_mochi_files = []
    for i in all_files:
        if i.endswith("_HOCD.mochi"):
            all_mochi_files.append(i[:-11])

    print("List of files:")
    for k in all_mochi_files:
        print(k)
    return all_mochi_files

#################################################################################################
def convert_ter_coords(input_file, output_file):
    """
    Convierte el archivo ter_coords.txt al formato requerido por el script.
    Utiliza la primera columna como ID del genoma y la tercera columna como coordenada de terminación.
    """
    print(f"Convirtiendo {input_file} al formato requerido...")
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            parts = line.strip().split(';')
            if len(parts) == 3:
                genome_id = parts[0].strip()
                terminus = parts[2].strip()  # Toma la tercera columna
                outfile.write(f"{genome_id}\t{terminus}\n")  # Formato requerido con tabulación
    print(f"Conversión completa. Archivo guardado como {output_file}.")

#################################################################################################
def load_ters(f):
    '''Carga las coordenadas de terminación del archivo convertido.'''
    crds = []
    with open(f, "r") as c:
        for line in c:
            splt = line.split("\t")
            crds.append((splt[0], int(splt[1].strip("\n"))))
    return crds

#################################################################################################
    
def match_ter(current_genome, all_ters):
    '''Empareja el genoma actual con su coordenada de terminación correspondiente.'''
    for entry in all_ters:
        # Permitir coincidencia parcial del nombre del genoma ignorando el sufijo
        if current_genome in entry[0]:
            return entry[1]
    return "Error"  # Loop should never end here

#################################################################################################
def get_coordinates(coords_file, ter):
    """Obtiene las coordenadas de los genes del archivo Mochiview."""
    coordinates_list = []
    linenum = 0
    mochilines = open(coords_file, "r").readlines()
    for line in mochilines:
        error = "OK"
        linenum += 1

        if linenum == 1:
            continue
        else:
            try:
                splt1 = line.split("\t")
                try:
                    start = int(splt1[1])
                    end = int(splt1[2])
                    if start > end:  # The first gene sometimes spans the origin.
                        start = 1
                except ValueError:  # Some mochiview lines will be "Error".
                    error = "Errant data point"

                ori = splt1[16].strip("\n")
                sys = splt1[4]
                strand = splt1[3]
                name = splt1[5]
                func = splt1[7]
                if ori == "HO":
                    if start < ter:
                        reverse_comp = "yes"
                    else:
                        reverse_comp = "no"
                else:
                    if start < ter:
                        reverse_comp = "no"
                    else:
                        reverse_comp = "yes"
            except IndexError:
                print("Error, mochi file is problematic. Likely that orientation is in the wrong column.")
                exit()
                
            if error == "OK":
                coordinates_list.append([start, end, ori, reverse_comp, name, sys, strand, func])
            else:
                continue

    return coordinates_list, mochilines

################################################################################################################################
def helper_RC(sequence):
    '''Returns the reverse-complement of the sequence.'''
    comp = ""
    rc = ""
    sequence_len = len(sequence) + 1

    for j in sequence:
        if j == "A" or j == "a":
            comp += "T"
        elif j == "C" or j == "c":
            comp += "G"
        elif j == "T" or j == "t":
            comp += "A"
        elif j == "G" or j == "g":
            comp += "C"
        else:
            comp += "A"
            print("Error during RC process: Bad Base -- ", j)
            print('Replacing bad base with "A"')
            continue

    for k in range(1, sequence_len):
        l = k * -1
        rc += comp[l]

    return rc

#################################################################################################
def analyze_sequence(crd, gnm_data):
    '''For each set of coordinates, gets nt sequence and counts the number of 3nt sequences.'''
    sequence_string = gnm_data[crd[0]-1:crd[1]]
    rc = crd[3]

    if rc == "yes":
        sequence_string = helper_RC(sequence_string)

    nts = [0, 0, 0, 0]  # number ACGT
    for nt in sequence_string:
        if nt == "A" or nt == "a":
            nts[0] += 1
        elif nt == "C" or nt == "c":
            nts[1] += 1
        elif nt == "G" or nt == "g":
            nts[2] += 1
        else:  # T
            nts[3] += 1

    return nts[0], nts[1], nts[2], nts[3], sequence_string

#################################################################################################
def print_analysis(cds, genomedata, original_mochi_data, genome_nm):
    '''Imprime los resultados en un archivo de salida.'''
    outfile = open(genome_nm + "_GC.mochi", "w")  # New mochi format file with GC skew column
    outfile.write(original_mochi_data[0].strip("\n") + "Orientation\t#A\t#T\t#C\t#G\tGCskew(rel. to LS)\tGCcontent\n" )

    forwig = []  # data will be returned as start, stop, gc_content
    HO = [0, 0, 0, 0]
    CD = [0, 0, 0, 0]

    place_counter = 0
    for entry in cds:
        place_counter += 1  # Mochiview line number index

        # Calculate GC skew
        ntA, ntC, ntG, ntT, gene_sequence = analyze_sequence(entry, genomedata)
        try:
            if entry[2] == "HO":
                skew = (ntC - ntG) / float(ntG + ntC)  # Antisense GC skew
            elif entry[2] == "CD":
                skew = (ntG - ntC) / float(ntG + ntC)  # GC skew
        except ZeroDivisionError:
            skew = 0.1

        # Record data
        GC_content = (ntG + ntC) / (ntG + ntA + ntC + ntT)
        forwig.append([entry[0], entry[1], skew, GC_content])
        new_mochi_line = original_mochi_data[place_counter].strip("\n") + "\t" + str(ntA) + "\t" + str(ntT) + "\t" + str(ntC) + "\t" + str(ntG) + "\t" + str(skew) + "\t" + str(GC_content) + "\n"
        outfile.write(new_mochi_line)
    outfile.close()

    return forwig

#################################################################################################
def genome_loader(g_file):
    '''Loads genome sequence.'''
    g = ""
    ctr = 0
    for l in open(g_file, "r"):
        ctr += 1
        if ctr > 1:
            g += l.strip("\n")
    return g

#################################################################################################
def GC_skew(seq, window):
    '''Calculates GC skew (G-C)/(G+C) for multiple windows along the sequence.'''
    GCcontentvalues = []
    values = []
    for i in range(0, len(seq), window):
        s = seq[i: i + window]
        g = s.count('G') + s.count('g')
        c = s.count('C') + s.count('c')
        a = s.count('A') + s.count('a')
        t = s.count('T') + s.count('t')
        try:
            skew = (g - c) / float(g + c)
            gcc = (g + c) / float(a + c + g + t)
        except ZeroDivisionError:
            skew = 0
            gcc = 0
        values.append(skew)
        GCcontentvalues.append(gcc)

    return values, GCcontentvalues

#################################################################################################
def collect_sequence(target):
    '''Collects just the sequence from a fasta file.'''
    seq = ""
    line_num = 0

    for line in open(target, "r"):
        line_num += 1

        if line_num == 1:
            continue
        else:
            seq += line.strip("\n")

    return seq

#################################################################################################
def skew_windows(vals, window):
    '''Exports GC values as a list of lists of "window-length step: gc value".'''
    data_steps = []
    line_ctr = 0

    for i in vals:
        line_ctr += 1
        x_point = line_ctr * window
        data_steps.append([x_point, (x_point + window - 1), i])

    return data_steps

#################################################################################################
def create_wig1(genome_len, gcdata, genomename, outfilename):
    '''Produces a new wig file for the whole genome.'''
    out = open(outfilename + ".wig", "w")
    out.write("track\nvariableStep chrom=" + genomename + " span=1\n")
    wigdata = []

    for i in range(1, genome_len):
        wigdata.append([i, 0])

    for j in gcdata:
        s = j[0]
        t = j[1]
        v = j[2]
        for k in range(s - 1, t):
            try:
                wigdata[k][1] = v
            except IndexError:
                continue

    for l in wigdata:
        out.write(str(l[0]) + "\t" + str(l[1]) + "\n")

    out.close()
    return None

#################################################################################################
def create_wig2(genome_len, gcdata, genomename, scale):
    '''Produces new wig files by gene for GC skew and GC content.'''
    out = open(genomename + "_GCS_gene.wig", "w")
    gcc_out = open(genomename + "_GCcont_gene.wig", "w")

    out.write("track\nvariableStep chrom=" + genomename + " span=1\n")
    gcc_out.write("track\nvariableStep chrom=" + genomename + " span=1\n")

    wigdata = []
    gccwig = []

    for i in range(1, genome_len):
        wigdata.append([i, 0])
        gccwig.append([i, 0])

    for j in gcdata:
        s = j[0]
        t = j[1]
        v = j[2]
        gc_cont = j[3]
        for k in range(s - 1, t):
            try:
                wigdata[k][1] = v
                gccwig[k][1] = gc_cont
            except IndexError:
                continue

    for l in wigdata:
        if l[1] != 0:
            out.write(str(l[0]) + "\t" + str(l[1] * scale) + "\n")
    for m in gccwig:
        if m[1] != 0:
            gcc_out.write(str(m[0]) + "\t" + str(m[1] * scale) + "\n")

    out.close()
    gcc_out.close()
    return None

#################################################################################################
# Function calls
#################################################################################################
print("Ensure that all files are in the same directory (for each genome):")
print("1. <genome_abbreviation>.fasta")
print("2. <genome_abbreviation>_HOCD.mochi")
print("Then execute like: $ python multi_allgeneGCskew.py")
answer = input("Continue? (Y/N): ")

if answer.lower() == "y":
    # Convertir el archivo ter_coords.txt al formato requerido
    convert_ter_coords("ter_coords.txt", "ter_coords_converted.txt")
    
    # Cargar el archivo convertido
    all_files = filenames()
    ter_data = load_ters("ter_coords_converted.txt")

    # Procesar cada genoma
    for a in all_files:
        print(f"Currently processing genome: {a}")
        mochi_file = f"{a}_HOCD.mochi"
        fastagenome = f"{a}.fasta"
        ter_position = match_ter(a, ter_data)
        if ter_position == "Error":
            print(f"Error: No terminus found for genome {a}")
            continue
        
        genomesequence = collect_sequence(fastagenome)
        genome_length = len(genomesequence) - 1

        # Calcula el GC skew para todo el genoma
        window_size = 100
        gcskew_values, GCCONT_VALS = GC_skew(genomesequence, window_size)
        data = skew_windows(gcskew_values, window_size)
        data2 = skew_windows(GCCONT_VALS, window_size)
        create_wig1(genome_length, data, a, (a + "_GCS_step"))
        create_wig1(genome_length, data2, a, (a + "_GCcont_step"))

        # Calcula el GC skew para genes individuales y escribe los datos en un nuevo archivo de Mochiview
        change_scale = 100
        cds1, mochiviewdata = get_coordinates(mochi_file, ter_position)
        data3 = print_analysis(cds1, genomesequence, mochiviewdata, a)
        create_wig2(genome_length, data3, a, change_scale)

    print("All processing has completed.")
else:
    print("Invalid entry. Exiting.")
    exit()


