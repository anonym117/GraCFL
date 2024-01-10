import matplotlib.pyplot as plt
import numpy as np
import os
import sys

def create_histo_scatter_xlog_log(reg_degree_list1, reg_degree_list2, bins_init, bins_final, graphname, file_path, analysis_type, freq_zero_init, freq_zero_final):
    degree_list_init = reg_degree_list1
    degree_list_final = reg_degree_list2
    
    plt.figure()
    xlabel0 = analysis_type + "-degree"
    plt.xlabel(xlabel0, fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.xscale('log')
    plt.yscale('log')

    # plt.xscale('log')
    # plt.yscale('log')
    # Create a figure and a set of subplots
    # fig, axs = plt.subplots(1, 2, figsize=(10, 4))  # 1 row, 2 columns

    title = analysis_type + "-degree distribution of " + graphname 
    # fig.suptitle(title, fontsize=16)
    plt.title(title, fontsize=14, pad=10)

    max_init = max(degree_list_init)
    min_init = min(degree_list_init)
    max_final = max(degree_list_final)
    min_final = min(degree_list_final)

    print (f"min_init: {min_init}")
    print (f"max_init: {max_init}")
    print (f"min_final: {min_final}")
    print (f"max_final: {max_final}")

    bin_edges1 = []
    bin_edges2 = []
    res = 1
    while res < max_init:
        bin_edges1.append(res)
        if res < 10:
            res += 1
        elif res < 100:
            res += 10
        elif res < 1000:
            res += 100
        elif res < 10000:
            res += 1000
        elif res < 100000:
            res += 10000
        elif res < 1000000:
            res += 100000
        elif res < 10000000:
            res += 1000000
    
    print(res)
    print(max_init)
    if res != max_init:
        bin_edges1.append(res)

    res = 1
    while res < max_final:
        bin_edges2.append(res)
        if res < 10:
            res += 1
        elif res < 100:
            res += 10
        elif res < 1000:
            res += 100
        elif res < 10000:
            res += 1000
        elif res < 100000:
            res += 10000
        elif res < 1000000:
            res += 100000
        elif res < 10000000:
            res += 1000000
    
    if res != max_final:
        bin_edges2.append(res)
    
    print(bin_edges1)
    print(bin_edges2)

    bin_edges1 = np.array(bin_edges1)
    bin_edges2 = np.array(bin_edges2)


    # bin_edges1 = np.arange(min(degree_list_init), max(degree_list_init) + 1, 1)
    # bin_edges2 = np.arange(min(degree_list_final), max(degree_list_final) + 1, 1)

    bin_counts1, _ = np.histogram(degree_list_init, bins=bin_edges1)
    bin_centers1 = 0.5 * (bin_edges1[:-1] + bin_edges1[1:])

    print(bin_counts1)
    print("bin_centers1")
    print(bin_centers1)

    print(f"bin_edge1 length: {len(bin_edges1)}")
    print(f"bin_centers1 length: {len(bin_centers1)}")
    print(f"bin_counts1 length: {len(bin_counts1)}")

    # bin_counts1 = np.insert(org_bin_counts1, 0, freq_zero_init)
    # bin_centers1 = np.insert(org_bin_centers1, 0, 0.1)

    # bin_counts1 = np.maximum(bin_counts1, 1e-5)

    bin_counts2, _ = np.histogram(degree_list_final, bins=bin_edges2)
    bin_centers2 = 0.5 * (bin_edges2[:-1] + bin_edges2[1:])
    # bin_counts2 = np.maximum(bin_counts2, 1e-5)

    print(bin_counts2)
    print("bin_centers2")
    print(bin_centers2)

    print(f"bin_edge2 length: {len(bin_edges2)}")
    print(f"bin_counts2 length: {len(bin_counts2)}")

    # bin_counts2 = np.insert(org_bin_counts2, 0, freq_zero_final)
    # bin_centers2 = np.insert(org_bin_centers2, 0, 0.12)

    # bin_centers1 = bin_edges1
    # bin_centers2 = bin_edges2

    plt.scatter(bin_centers1, bin_counts1, color='blue', marker='s', label='Initial Graph', s=10)
    plt.scatter(bin_centers2, bin_counts2, color='green', marker='x', label='Final Graph', s=10)

    xmin, xmax = plt.xlim()

    # print(bin_counts1)
    # print(bin_counts2)
    # print(bin_centers1)
    # print(bin_centers2)

    # Adjust layout for better visualization
    plt.tight_layout()

    # # fig.subplots_adjust(top=0.85)
    # # Add a legend
    plt.legend()

    # # Show the plot
    plt.show()

    output_filename_pdf = "histogram-" + graphname + "-" + analysis_type + "-dot-xlog-1-histo-all-5-log.pdf"
    output_filename_png = "histogram-" + graphname + "-" + analysis_type + "-dot-xlog-1-histo-all-5-log.png"

    if (os.path.exists(output_filename_pdf)):
        os.remove(output_filename_pdf)

    if (os.path.exists(output_filename_png)):
        os.remove(output_filename_png)
    plt.savefig(output_filename_pdf)
    plt.savefig(output_filename_png)

def create_histo_scatter_xlog_all(reg_degree_list1, reg_degree_list2, bins_init, bins_final, graphname, file_path, analysis_type, freq_zero_init, freq_zero_final):
    # print(type(degree_list))
    # degree_list = np.array(reg_degree_list)
    degree_list_init = reg_degree_list1
    degree_list_final = reg_degree_list2
    
    plt.figure()
    xlabel0 = analysis_type + "-degree"
    plt.xlabel(xlabel0, fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.xscale('log')
    plt.yscale('log')

    # plt.xscale('log')
    # plt.yscale('log')
    # Create a figure and a set of subplots
    # fig, axs = plt.subplots(1, 2, figsize=(10, 4))  # 1 row, 2 columns

    title = analysis_type + "-degree distribution of " + graphname 
    # fig.suptitle(title, fontsize=16)
    plt.title(title, fontsize=14, pad=10)

    bin_edges1 = np.arange(1, max(degree_list_init) + 1, 1)
    bin_edges2 = np.arange(1, max(degree_list_final) + 1, 1)

    print(f"bin_edges1 {bin_edges1}")
    print(f"bin_edges2 {bin_edges2}")

    bin_counts1, _ = np.histogram(degree_list_init, bins=bin_edges1)
    bin_centers1 = 0.5 * (bin_edges1[:-1] + bin_edges1[1:])

    # print(org_bin_counts1)
    # print(org_bin_centers1)

    # bin_counts1 = np.insert(org_bin_counts1, 0, freq_zero_init)
    # bin_centers1 = np.insert(org_bin_centers1, 0, 0.1)

    # bin_counts1 = np.maximum(bin_counts1, 1e-5)

    bin_counts2, _ = np.histogram(degree_list_final, bins=bin_edges2)
    bin_centers2 = 0.5 * (bin_edges2[:-1] + bin_edges2[1:])
    # bin_counts2 = np.maximum(bin_counts2, 1e-5)

    # print(org_bin_counts2)
    # print(org_bin_centers2)

    # bin_counts2 = np.insert(org_bin_counts2, 0, freq_zero_final)
    # bin_centers2 = np.insert(org_bin_centers2, 0, 0.12)

    plt.scatter(bin_centers1, bin_counts1, color='blue', marker='s', label='Initial Graph', s=5)
    plt.scatter(bin_centers2, bin_counts2, color='green', marker='x', label='Final Graph', s=5)

    xmin, xmax = plt.xlim()

    # print(bin_counts1)
    # print(bin_counts2)
    # print(bin_centers1)
    # print(bin_centers2)sss

    # Adjust layout for better visualization
    plt.tight_layout()

    # fig.subplots_adjust(top=0.85)
    # Add a legend
    plt.legend()

    # Show the plot
    plt.show()

    output_filename_pdf = "histogram-" + graphname + "-" + analysis_type + "-dot-xlog-1-histo-all-5.pdf"
    output_filename_png = "histogram-" + graphname + "-" + analysis_type + "-dot-xlog-1-histo-all-5.png"

    if (os.path.exists(output_filename_pdf)):
        os.remove(output_filename_pdf)

    if (os.path.exists(output_filename_png)):
        os.remove(output_filename_png)
    plt.savefig(output_filename_pdf)
    plt.savefig(output_filename_png)

def create_histo_scatter_xlog_2_histo(reg_degree_list1, reg_degree_list2, bins_init, bins_final, graphname, file_path, analysis_type, freq_zero_init, freq_zero_final):
    # print(type(degree_list))
    # degree_list = np.array(reg_degree_list)
    degree_list_init = reg_degree_list1
    degree_list_final = reg_degree_list2
    
    plt.figure()
    xlabel0 = analysis_type + "-degree"
    plt.xlabel(xlabel0, fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.xscale('log')
    plt.yscale('log')

    # plt.xscale('log')
    # plt.yscale('log')
    # Create a figure and a set of subplots
    # fig, axs = plt.subplots(1, 2, figsize=(10, 4))  # 1 row, 2 columns

    title = analysis_type + "-degree distribution of " + graphname 
    # fig.suptitle(title, fontsize=16)
    plt.title(title, fontsize=14, pad=10)

    org_bin_counts1, bin_edges1 = np.histogram(degree_list_init, bins=bins_init)
    org_bin_centers1 = 0.5 * (bin_edges1[:-1] + bin_edges1[1:])

    # print(org_bin_counts1)
    # print(org_bin_centers1)

    bin_counts1 = np.insert(org_bin_counts1, 0, freq_zero_init)
    bin_centers1 = np.insert(org_bin_centers1, 0, 0.1)

    # bin_counts1 = np.maximum(bin_counts1, 1e-5)

    org_bin_counts2, bin_edges2 = np.histogram(degree_list_final, bins=bins_final)
    org_bin_centers2 = 0.5 * (bin_edges2[:-1] + bin_edges2[1:])
    # bin_counts2 = np.maximum(bin_counts2, 1e-5)

    # print(org_bin_counts2)
    # print(org_bin_centers2)

    bin_counts2 = np.insert(org_bin_counts2, 0, freq_zero_final)
    bin_centers2 = np.insert(org_bin_centers2, 0, 0.12)

    plt.scatter(bin_centers1, bin_counts1, color='blue', label='Initial Graph')
    plt.scatter(bin_centers2, bin_counts2, color='green', label='Final Graph')

    xmin, xmax = plt.xlim()

    # print(bin_counts1)
    # print(bin_counts2)
    # print(bin_centers1)
    # print(bin_centers2)

    # Adjust the minimum x-axis limit to reduce space
    # plt.xlim(min(bin_counts1), xmax)

    # Plot histogram on the first subplot
    # plt.scatter(bin_centers, bin_counts, color='blue')
    # axs[0].scatter(bin_centers1, bin_counts1, color='blue')
    # # axs[0].hist(data1, bins=30, color='blue', alpha=0.7)
    # axs[0].set_title('Initial Graph')
    # xlabel0 = analysis_type + "-degree"
    # axs[0].set_xlabel(xlabel0)
    # axs[0].set_ylabel('Frequency')
    # axs[0].set_xscale('log')
    # axs[0].set_yscale('log')
    # axs[0].xscale('log')
    # axs[0].yscale('log')

    # Plot histogram on the second subplot
    # axs[1].scatter(bin_centers2, bin_counts2, color='green')
    # # axs[1].hist(data2, bins=30, color='green', alpha=0.7)
    # axs[1].set_title('Final Graph')
    # xlabel1 = analysis_type + "-degree"
    # axs[1].set_xlabel(xlabel1)
    # axs[1].set_ylabel('Frequency')
    # axs[1].set_xscale('log')
    # axs[1].set_yscale('log')

    # Adjust layout for better visualization
    plt.tight_layout()

    # fig.subplots_adjust(top=0.85)
    # Add a legend
    plt.legend()

    # Show the plot
    plt.show()

    output_filename_pdf = "histogram-" + graphname + "-" + analysis_type + "-dot-xlog-1-histo.pdf"
    output_filename_png = "histogram-" + graphname + "-" + analysis_type + "-dot-xlog-1-histo.png"

    if (os.path.exists(output_filename_pdf)):
        os.remove(output_filename_pdf)

    if (os.path.exists(output_filename_png)):
        os.remove(output_filename_png)
    plt.savefig(output_filename_pdf)
    plt.savefig(output_filename_png)


    # plt.savefig('test-4-2histo.pdf')

    # plt.hist(degree_list, bins=bins, alpha=0.85, color='blue', edgecolor='black')
    # plt.yscale('log')
    # plt.xscale('log')

    # Add titles and labels
    # plt.title(title)
    # xlabel = analysis_type+"-degree"
    # plt.xlabel(xlabel)
    # plt.ylabel('Frequency')
    # # plt.tight_layout()

    # # Get the current x-axis limits
    # xmin, xmax = plt.xlim()

    # # Adjust the minimum x-axis limit to reduce space
    # plt.xlim(min(degree_list), max(degree_list)+1)

    # bin_counts, bin_edges = np.histogram(degree_list, bins=bins)
    # # print(bin_edges)
    # # print(bin_counts)
    # bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    # # Plot each bin's frequency as a dot
    # plt.scatter(bin_centers, bin_counts, color='blue')

    # output_filename = "histogram-" + file_path + "-dot-xlog.pdf"
    # if (os.path.exists(output_filename)):
    #     os.remove(output_filename)

    # plt.savefig(output_filename)

    # # Show the plot
    # plt.show()

    # fig, axs = plt.subplots(1, 2, figsize=(10, 4))  # 1 row, 2 columns
    # fig.suptitle('Comparison of Two Datasets', fontsize=16)

    # # Plot histogram on the first subplot
    # axs[0].hist(degree_list1, bins=30, color='blue', alpha=0.7)
    # axs[0].set_title('Histogram of Data 1')
    # axs[0].set_xlabel('Value')
    # axs[0].set_ylabel('Frequency')

    # # Plot histogram on the second subplot
    # axs[1].hist(reg_degree_list2, bins=30, color='green', alpha=0.7)
    # axs[1].set_title('Histogram of Data 2')
    # axs[1].set_xlabel('Value')
    # axs[1].set_ylabel('Frequency')

    # # Adjust layout for better visualization
    # plt.tight_layout()

    # # Show the plot
    # plt.show()

def process_degree_to_vtx(lines_init, file_path_init, lines_final, file_path_final, analysis_type, graphname):
    degree_list_init = []
    degree_list_final = []

    freq_zero_init = 0

    for line in lines_init:
        line_splits = line.strip().split(" ")
        for degree in line_splits:
            try:
                degree_list_init.append(int(degree.strip()))
                if int(degree.strip()) >0:
                    degree_list_init.append(int(degree.strip()))
                else:
                    freq_zero_init = freq_zero_init + 1
                # if int(degree.strip()) == 0:
                #     freq_zero_init += 1
            except ValueError as e:
                print(f"Could not convert '{degree}' to an integer: {e}")

    # for line in lines_init:
    #     line_splits = line.strip().split(" ")
    #     for degree in line_splits:
    #         try:
    #             degree_list_init.append(int(degree.strip()))
    #         except ValueError as e:
    #             print(f"Could not convert '{degree}' to an integer: {e}")

    # print(degree_list)
    # Create the histogram
    max_degree_init = max(degree_list_init)
    bins_init = 100
    if (max_degree_init < 100):
        bins_init = max_degree_init

    freq_zero_final = 0

    for line in lines_final:
        line_splits = line.strip().split(" ")
        for degree in line_splits:
            try:
                degree_list_final.append(int(degree.strip()))
                if int(degree.strip()) >0:
                    degree_list_final.append(int(degree.strip()))
                else:
                    freq_zero_final = freq_zero_final + 1
                # if int(degree.strip()) == 0:
                #     freq_zero_final += 1
            except ValueError as e:
                print(f"Could not convert '{degree}' to an integer: {e}")

    # print(degree_list)
    # Create the histogram
    max_degree_final = max(degree_list_final)
    bins_final = 100
    if (max_degree_final < 100):
        bins_final = max_degree_final

    print(f"freq_zero_init: {freq_zero_init}")
    print(f"freq_zero_final: {freq_zero_final}")

    create_histo_scatter_xlog_log(degree_list_init, degree_list_final, bins_init, bins_final, graphname, file_path_init, analysis_type, freq_zero_init, freq_zero_final)
    # create_histo_scatter_xlog_2_histo(degree_list_init, degree_list_final, bins_init, bins_final, graphname, file_path_init, analysis_type, freq_zero_init, freq_zero_final)
    create_histo_scatter_xlog_all(degree_list_init, degree_list_final, bins_init, bins_final, graphname, file_path_init, analysis_type, freq_zero_init, freq_zero_final)
    # create_histo(degree_list, bins, graphname, file_path, analysis_type)
    # create_histo_scatter(degree_list, bins, graphname, file_path, analysis_type)
    # create_histo_scatter_unequal(degree_list, bins, graphname, file_path, analysis_type)
    

def main():
    # Replace 'yourfile.txt' with the path to the file you want to read
    file_path_init = 'vToDeg-out-httpd_df-init'
    file_path_final = 'vToDeg-out-httpd_df-final'
    graphname = ""
    analysis_type = ""
    if len(sys.argv) == 5:
        file_path_init = str(sys.argv[1])
        file_path_final = str(sys.argv[2])
        # ex: httpd_df
        graphname = str(sys.argv[3])
        # analysis_type: in/out
        analysis_type = str(sys.argv[4]) 
    print (file_path_init)
    print (file_path_final)

    try:
        file_init = open(file_path_init, 'r')
        lines_init = file_init.readlines()
        file_init.close()

        file_final = open(file_path_final, 'r')
        lines_final = file_final.readlines()
        file_final.close()
        # process_degree_to_vtx(lines, file_path, analysis_type, graphname)
        process_degree_to_vtx(lines_init, file_path_init, lines_final, file_path_final, analysis_type, graphname)
        # process_degree_to_vtx_without_zero(lines, file_path, analysis_type, graphname)
            
    except FileNotFoundError:
        print("File not found.")
    except IOError:
        print("Error while reading the file.")

if __name__ == "__main__":
    main()