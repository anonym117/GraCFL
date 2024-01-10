import matplotlib.pyplot as plt
import numpy as np
import os
import sys

def create_histo_scatter_xlog_2_histo(reg_degree_list1, reg_degree_list2, bins, graphname, file_path, analysis_type):
    # print(type(degree_list))
    # degree_list = np.array(reg_degree_list)
    degree_list_init = reg_degree_list1
    degree_list_final = reg_degree_list2
    
    # plt.figure()
    # Create a figure and a set of subplots
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))  # 1 row, 2 columns

    plt.yscale('log')
    plt.xscale('log')

    title = "Histogram - " + analysis_type + " degree distribution of " + graphname 
    fig.suptitle(title, fontsize=16)

    bin_counts1, bin_edges1 = np.histogram(degree_list_init, bins=100)
    bin_centers1 = 0.5 * (bin_edges1[:-1] + bin_edges1[1:])

    bin_counts2, bin_edges2 = np.histogram(degree_list_final, bins=100)
    bin_centers2 = 0.5 * (bin_edges2[:-1] + bin_edges2[1:])

    # Plot histogram on the first subplot
    # plt.scatter(bin_centers, bin_counts, color='blue')
    axs[0].scatter(bin_centers1, bin_counts1, color='blue')
    # axs[0].hist(data1, bins=30, color='blue', alpha=0.7)
    axs[0].set_title('Histogram of Data 1')
    axs[0].set_xlabel('Value')
    axs[0].set_ylabel('Frequency')

    # Plot histogram on the second subplot
    axs[1].scatter(bin_centers2, bin_counts2, color='green')
    # axs[1].hist(data2, bins=30, color='green', alpha=0.7)
    axs[1].set_title('Histogram of Data 2')
    axs[1].set_xlabel('Value')
    axs[1].set_ylabel('Frequency')

    # Adjust layout for better visualization
    plt.tight_layout()

    # Show the plot
    plt.show()

    output_filename_pdf = "histogram-" + file_path + "-dot-xlog-2-histo.pdf"
    output_filename_png = "histogram-" + file_path + "-dot-xlog-2-histo.png"

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

def create_histo_scatter_xlog(reg_degree_list, bins, graphname, file_path, analysis_type):
    # print(type(degree_list))
    # degree_list = np.array(reg_degree_list)
    degree_list = reg_degree_list
    
    plt.figure()
    # Create a figure and a set of subplots
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))  # 1 row, 2 columns

    title = "Histogram - " + analysis_type + " degree distribution of " + graphname 
    fig.suptitle('Comparison of Two Datasets', fontsize=16)
    # plt.hist(degree_list, bins=bins, alpha=0.85, color='blue', edgecolor='black')
    plt.yscale('log')
    plt.xscale('log')

    # Add titles and labels
    plt.title(title)
    xlabel = analysis_type+"-degree"
    plt.xlabel(xlabel)
    plt.ylabel('Frequency')
    # plt.tight_layout()

    # Get the current x-axis limits
    xmin, xmax = plt.xlim()

    # Adjust the minimum x-axis limit to reduce space
    plt.xlim(min(degree_list), max(degree_list)+1)

    bin_counts, bin_edges = np.histogram(degree_list, bins=bins)
    # print(bin_edges)
    # print(bin_counts)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    # Plot each bin's frequency as a dot
    plt.scatter(bin_centers, bin_counts, color='blue')

    output_filename = "histogram-" + file_path + "-dot-xlog.pdf"
    if (os.path.exists(output_filename)):
        os.remove(output_filename)

    plt.savefig(output_filename)

    # Show the plot
    plt.show()


def create_histo_scatter_unequal(reg_degree_list, bins, graphname, file_path, analysis_type):
    plt.figure()
    # print(reg_degree_list)
    # degree_list = np.array(reg_degree_list)
    degree_list = reg_degree_list
    # print(degree_list)
    # plt.hist(degree_list, bins=bins, alpha=0.85, color='blue', edgecolor='black')
    min_elem = min(degree_list)
    max_elem = max(degree_list)

    start_elem = 0
    is_less = 0
    end_elem = max(degree_list)
    bin_edges = []

    if end_elem < 100:
        is_less = -1

    if is_less == -1:
        for i in range(0, end_elem+1, 1):
            bin_edges.append(i)
    else:
        for i in range(0, 20, 1):
            bin_edges.append(i)

        total_range = end_elem + 1 - 20
        bin_eq_range = 0.0125 * total_range
        # print(bin_eq_range)
        # j = -1
        # for i in range(20, end_elem + 1, bin_eq_range):
        #     bin_edges.append(i)
        i = 20
        while i < end_elem+1.0:
            bin_edges.append(int(i))
            i += bin_eq_range
    
    # print (bin_edges)
    # print(end_elem)

    plt.yscale('log')
    # plt.xscale('log')

    # Add titles and labels
    title = "Histogram - " + analysis_type + " degree distribution of " + graphname 
    plt.title(title)
    xlabel = analysis_type+"-degree"
    plt.xlabel(xlabel)
    plt.ylabel('Frequency')
    # plt.tight_layout()

    # Get the current x-axis limits
    xmin, xmax = plt.xlim()

    # Adjust the minimum x-axis limit to reduce space
    # plt.xlim(min(degree_list), xmax)
    plt.xlim(min(degree_list)-1, max(degree_list)+1)

    bin_counts, bin_edges = np.histogram(degree_list, bins=bin_edges)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    # print(bin_centers)

    # Plot each bin's frequency as a dot
    plt.scatter(bin_centers, bin_counts, color='blue')

    output_filename = "histogram-" + file_path + "-dot-unequal.pdf"
    if (os.path.exists(output_filename)):
        os.remove(output_filename)

    plt.savefig(output_filename)
    print(output_filename)

    # Show the plot
    plt.show()


def create_histo_scatter_unequal_xlog(reg_degree_list, bins, graphname, file_path, analysis_type):
    plt.figure()
    # print(reg_degree_list)
    # degree_list = np.array(reg_degree_list)
    degree_list = reg_degree_list
    # print(degree_list)
    # plt.hist(degree_list, bins=bins, alpha=0.85, color='blue', edgecolor='black')
    min_elem = min(degree_list)
    max_elem = max(degree_list)

    start_elem = 0
    is_less = 0
    end_elem = max(degree_list)
    bin_edges = []

    if end_elem < 100:
        is_less = -1

    if is_less == -1:
        for i in range(0, end_elem+1, 1):
            bin_edges.append(i)
    else:
        for i in range(0, 20, 1):
            bin_edges.append(i)

        total_range = end_elem + 1 - 20
        bin_eq_range = 0.0125 * total_range
        print(bin_eq_range)
        # j = -1
        # for i in range(20, end_elem + 1, bin_eq_range):
        #     bin_edges.append(i)
        i = 20
        while i < end_elem+1.0:
            bin_edges.append(int(i))
            i += bin_eq_range
    
    print (bin_edges)
    print(end_elem)

    plt.yscale('log')
    plt.xscale('log')

    # Add titles and labels
    title = "Histogram - " + analysis_type + " degree distribution of " + graphname 
    plt.title(title)
    xlabel = analysis_type+"-degree"
    plt.xlabel(xlabel)
    plt.ylabel('Frequency')
    # plt.tight_layout()

    # Get the current x-axis limits
    xmin, xmax = plt.xlim()

    # Adjust the minimum x-axis limit to reduce space
    # plt.xlim(min(degree_list), xmax)
    plt.xlim(min(degree_list)-1, max(degree_list)+1)

    bin_counts, bin_edges = np.histogram(degree_list, bins=bin_edges)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    # print(bin_centers)

    # Plot each bin's frequency as a dot
    plt.scatter(bin_centers, bin_counts, color='blue')

    output_filename = "histogram-" + file_path + "-dot-unequal-xlog.pdf"
    if (os.path.exists(output_filename)):
        os.remove(output_filename)

    plt.savefig(output_filename)
    print(output_filename)

    # Show the plot
    plt.show()


def create_histo_scatter(reg_degree_list, bins, graphname, file_path, analysis_type):
    plt.figure()
    # print(reg_degree_list)
    # degree_list = np.array(reg_degree_list)
    degree_list = reg_degree_list
    # print(degree_list)
    # plt.hist(degree_list, bins=bins, alpha=0.85, color='blue', edgecolor='black')
    plt.yscale('log')
    # plt.xscale('log')

    # Add titles and labels
    title = "Histogram - " + analysis_type + " degree distribution of " + graphname 
    plt.title(title)
    xlabel = analysis_type+"-degree"
    plt.xlabel(xlabel)
    plt.ylabel('Frequency')
    # plt.tight_layout()

    # Get the current x-axis limits
    xmin, xmax = plt.xlim()

    # Adjust the minimum x-axis limit to reduce space
    # plt.xlim(min(degree_list), xmax)
    plt.xlim(min(degree_list)-1, max(degree_list)+1)

    bin_counts, bin_edges = np.histogram(degree_list, bins=bins)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    # print(bin_counts)
    # print(bin_centers)

    # print(bin_centers)

    # Plot each bin's frequency as a dot
    plt.scatter(bin_centers, bin_counts, color='blue')

    output_filename = "histogram-" + file_path + "-dot.pdf"
    if (os.path.exists(output_filename)):
        os.remove(output_filename)

    plt.savefig(output_filename)
    print(output_filename)

    # Show the plot
    plt.show()


def create_histo_scatter_xlog(reg_degree_list, bins, graphname, file_path, analysis_type):
    # print(type(degree_list))
    # degree_list = np.array(reg_degree_list)
    degree_list = reg_degree_list
    
    plt.figure()
    # plt.hist(degree_list, bins=bins, alpha=0.85, color='blue', edgecolor='black')
    plt.yscale('log')
    plt.xscale('log')

    # Add titles and labels
    title = "Histogram - " + analysis_type + " degree distribution of " + graphname 
    plt.title(title)
    xlabel = analysis_type+"-degree"
    plt.xlabel(xlabel)
    plt.ylabel('Frequency')
    # plt.tight_layout()

    # Get the current x-axis limits
    xmin, xmax = plt.xlim()

    # Adjust the minimum x-axis limit to reduce space
    plt.xlim(min(degree_list), max(degree_list)+1)

    bin_counts, bin_edges = np.histogram(degree_list, bins=bins)
    # print(bin_edges)
    # print(bin_counts)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    # Plot each bin's frequency as a dot
    plt.scatter(bin_centers, bin_counts, color='blue')

    output_filename = "histogram-" + file_path + "-dot-xlog.pdf"
    if (os.path.exists(output_filename)):
        os.remove(output_filename)

    plt.savefig(output_filename)

    # Show the plot
    plt.show()

def create_histo(degree_list, bins, graphname, file_path, analysis_type):
    plt.figure()
    plt.hist(degree_list, bins=bins, alpha=0.85, color='blue', edgecolor='black')
    plt.yscale('log')
    # plt.xscale('log')

    # Add titles and labels
    title = "Histogram - " + analysis_type + " degree distribution of " + graphname 
    plt.title(title)
    xlabel = analysis_type+"-degree"
    plt.xlabel(xlabel)
    plt.ylabel('Frequency')
    plt.tight_layout()

    # Get the current x-axis limits
    xmin, xmax = plt.xlim()

    # Adjust the minimum x-axis limit to reduce space
    plt.xlim(min(degree_list), xmax)

    output_filename = "histogram-" + file_path + ".pdf"
    if (os.path.exists(output_filename)):
        os.remove(output_filename)

    plt.savefig(output_filename)

    # Show the plot
    plt.show()

def create_histo_xlog(degree_list, bins, graphname, file_path, analysis_type):
    plt.figure()
    plt.hist(degree_list, bins=bins, alpha=0.85, color='blue', edgecolor='black', histtype='step', linewidth=2)
    plt.yscale('log')
    plt.xscale('log')
    # plt.xscale('log')

    # Add titles and labels
    title = "Histogram - " + analysis_type + " degree distribution of " + graphname 
    plt.title(title)
    xlabel = analysis_type+"-degree"
    plt.xlabel(xlabel)
    plt.ylabel('Frequency')
    plt.tight_layout()

    # Get the current x-axis limits
    xmin, xmax = plt.xlim()

    # Adjust the minimum x-axis limit to reduce space
    plt.xlim(min(degree_list), xmax)

    output_filename = "histogram-" + file_path + "-without_zero.pdf"
    if (os.path.exists(output_filename)):
        os.remove(output_filename)

    plt.savefig(output_filename)

    # Show the plot
    plt.show()



def process_degree_to_vtx(lines, file_path, analysis_type, graphname):
    degree_list = []

    for line in lines:
        line_splits = line.strip().split(" ")
        for degree in line_splits:
            try:
                degree_list.append(int(degree.strip()))
            except ValueError as e:
                print(f"Could not convert '{degree}' to an integer: {e}")

    # print(degree_list)
    # Create the histogram
    max_degree = max(degree_list)
    bins = 100
    if (max_degree < 100):
        bins = max_degree

    create_histo(degree_list, bins, graphname, file_path, analysis_type)
    create_histo_scatter(degree_list, bins, graphname, file_path, analysis_type)
    create_histo_scatter_unequal(degree_list, bins, graphname, file_path, analysis_type)


def process_degree_to_vtx_without_zero(lines, file_path, analysis_type, graphname):
    degree_list_without_zero = []
    freq_zero = 0

    for line in lines:
        line_splits = line.strip().split(" ")
        for degree in line_splits:
            try:
                if int(degree.strip()) >0:
                    degree_list_without_zero.append(int(degree.strip()))
                else:
                    freq_zero += 1
            except ValueError as e:
                print(f"Could not convert '{degree}' to an integer: {e}")

    # print(degree_list)
    # Create the histogram
    max_degree = max(degree_list_without_zero)
    bins = 100
    if (max_degree < 100):
        bins = max_degree

    create_histo_xlog(degree_list_without_zero, bins, graphname, file_path, analysis_type)
    create_histo_scatter_xlog(degree_list_without_zero, bins, graphname, file_path, analysis_type)
    create_histo_scatter_unequal_xlog(degree_list_without_zero, bins, graphname, file_path, analysis_type)

    print("Zero frequency")
    print(freq_zero)
    


    

def main():
    # Replace 'yourfile.txt' with the path to the file you want to read
    file_path = 'vToDeg-out-httpd_df-final'
    graphname = ""
    analysis_type = ""
    if len(sys.argv) == 4:
        file_path = str(sys.argv[1])
        # ex: httpd_df
        graphname = str(sys.argv[2])
        # analysis_type: in/out
        analysis_type = str(sys.argv[3]) 
    print (file_path)

    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
            process_degree_to_vtx(lines, file_path, analysis_type, graphname)
            process_degree_to_vtx_without_zero(lines, file_path, analysis_type, graphname)
            
    except FileNotFoundError:
        print("File not found.")
    except IOError:
        print("Error while reading the file.")

if __name__ == "__main__":
    main()