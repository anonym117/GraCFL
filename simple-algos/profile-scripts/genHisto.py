import matplotlib.pyplot as plt
import os


def generate_histo(degree_to_vtx_dict_init, max_degree_init, title,  xlabel, ylabel, bars_no, filename):
    histo_range_x_init = max_degree_init // bars_no
    # print(histo_range_x_init)
    categories_dic_init = {}
    bin_edges_init = []
    start_init = 0


    finish_before = max_degree_init + 1 - histo_range_x_init
    
    for i in range(1, finish_before, histo_range_x_init):
        key = "" + str(start_init) + "-" + str(start_init + histo_range_x_init) + ""

        value = 0
        for j in range(start_init, start_init+ histo_range_x_init):
            if str(j) in degree_to_vtx_dict_init:
                value += int(degree_to_vtx_dict_init[str(j)])

        bin_edges_init.append(start_init)
        start_init += histo_range_x_init

        categories_dic_init[key] = value

    key = "" + str(finish_before) + "-" + str(max_degree_init) +""
    value = 0
    for j in range(finish_before, max_degree_init + 1):
        if str(j) in degree_to_vtx_dict_init:
            value += int(degree_to_vtx_dict_init[str(j)])
    categories_dic_init[key] = value
    #bin_edges_init.append(max_degree_final)
    



    # Names of categories and their frequencies
    categories_init = list(categories_dic_init.keys())
    frequencies_init = list(categories_dic_init.values())

    print("------------")
    print(max_degree_init)
    print(categories_init)
    print(frequencies_init)

    plt.figure()
    # Making the y-axis logarithmic
    plt.yscale('log')
    # Creating the histogram
    # plt.hist(frequencies_init, bins=bin_edges_init, edgecolor='black')
    #This is another way of creating histogram
    bars = plt.bar(categories_init, frequencies_init)

    for bar in bars:
        height = bar.get_height()
        formatted_height = "{:.1e}".format(height)
        height_value = float(formatted_height.split('e')[0])
        print(height_value)
        plt.annotate(f'{height_value}', 
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3),  # 3 points vertical offset
                textcoords="offset points",
                ha='center', va='bottom', fontsize=6)

    # bin_labels = categories_init
    # plt.xticks([(a + b) / 2 for a, b in zip(bin_edges_init[:-1], bin_edges_init[1:])], bin_labels)

    # Adding titles and labels
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xticks(rotation=45, fontsize=10)
    plt.yticks(fontsize=10)
    # Adjust layout and bottom margin
    plt.tight_layout(pad=2)
    plt.subplots_adjust(bottom=0.25)  # Adjust bottom margin
    # plt.subplots_adjust(top=0.85)


    # Save the plot to a file
    plt.savefig(filename)  # Saves the histogram as a PNG file
    # plt.savefig('histogram-degree-to-vtx-init.pdf')  # Uncomment to save as a PDF

    # Display the histogram
    plt.show()


def generate_histo_2(degree_to_vtx_dict_init, max_degree_init, title,  xlabel, ylabel, bars_no, filename):
    histo_range_x_init = max_degree_init // bars_no
    # print(histo_range_x_init)
    categories_dic_init = {}
    bin_edges_init = []
    start_init = 0


    finish_before = max_degree_init + 1 - histo_range_x_init
    
    for i in range(1, finish_before, histo_range_x_init):
        key = start_init + histo_range_x_init

        value = 0
        for j in range(start_init, start_init+ histo_range_x_init):
            if str(j) in degree_to_vtx_dict_init:
                value += int(degree_to_vtx_dict_init[str(j)])

        bin_edges_init.append(start_init + histo_range_x_init)
        categories_dic_init[key] = value

        start_init += histo_range_x_init

    key =  max_degree_init
    value = 0
    for j in range(finish_before, max_degree_init + 1):
        if str(j) in degree_to_vtx_dict_init:
            value += int(degree_to_vtx_dict_init[str(j)])
    categories_dic_init[key] = value
    bin_edges_init.append(max_degree_init)
    
    # Names of categories and their frequencies
    categories_init = list(categories_dic_init.keys())
    frequencies_init = list(categories_dic_init.values())

    print("------------")
    print(max_degree_init)
    print(categories_init)
    print(frequencies_init)

    plt.figure()
    # Making the y-axis logarithmic
    plt.yscale('log')
    # Creating the histogram
    # plt.hist(frequencies_init, bins=bin_edges_init, edgecolor='black')
    #This is another way of creating histogram
    bars = plt.bar(categories_init, frequencies_init)

    for bar in bars:
        height = bar.get_height()
        formatted_height = "{:.1e}".format(height)
        height_value = float(formatted_height.split('e')[0])
        print(height_value)
        plt.annotate(f'{height_value}', 
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3),  # 3 points vertical offset
                textcoords="offset points",
                ha='center', va='bottom', fontsize=6)

    # bin_labels = categories_init
    # plt.xticks([(a + b) / 2 for a, b in zip(bin_edges_init[:-1], bin_edges_init[1:])], bin_labels)

    # Adding titles and labels
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xticks(rotation=45, fontsize=10)
    plt.yticks(fontsize=10)
    # Adjust layout and bottom margin
    plt.tight_layout(pad=2)
    plt.subplots_adjust(bottom=0.25)  # Adjust bottom margin
    # plt.subplots_adjust(top=0.85)


    # Save the plot to a file
    plt.savefig(filename)  # Saves the histogram as a PNG file
    # plt.savefig('histogram-degree-to-vtx-init.pdf')  # Uncomment to save as a PDF

    # Display the histogram
    plt.show()


def process_out_degree_to_vtx(lines):
    # outDegree = 0, vtxCnt = 417565
    flag = 0
    degree_to_vtx_dict_init = {}
    degree_to_vtx_dict_final = {}
    max_degree_init = -1
    max_degree_final = -1
    
    for  i in range(len(lines)):
        if lines[i] == "----outDegreeToVtx----\n" and flag == 0:
            flag = 1
            j = i+1
            while lines[j] != "----inDegreeToVtx----\n":
                line_splits = lines[j].strip().split(",")
                # print(line_splits)
                out_degree = line_splits[0].split(" = ")[1].strip()
                if max_degree_init < int(out_degree):
                    max_degree_init = int(out_degree)
                vtx_cnt = line_splits[1].split(" = ")[1].strip()
                if out_degree in degree_to_vtx_dict_init:
                    degree_to_vtx_dict_init[out_degree] += vtx_cnt
                else:
                    degree_to_vtx_dict_init[out_degree] = vtx_cnt
                j = j+1
                # print(j)
        
        if lines[i] == "----outDegreeToVtx----\n" and flag == 1:
            j = i+1
            while lines[j] != "----inDegreeToVtx----\n":
                line_splits = lines[j].strip().split(",")
                # print(line_splits)
                out_degree = line_splits[0].split(" = ")[1].strip()
                if max_degree_final < int(out_degree):
                    max_degree_final = int(out_degree)
                vtx_cnt = line_splits[1].split(" = ")[1].strip()
                if out_degree in degree_to_vtx_dict_final:
                    degree_to_vtx_dict_final[out_degree] += vtx_cnt
                else:
                    degree_to_vtx_dict_final[out_degree] = vtx_cnt
                j = j+1
                # print(j)

    # for key, value in degree_to_vtx_dict_init.items():
    #     print(key + ", " + value + "\n")
    
    # for key, value in degree_to_vtx_dict_final.items():
    #     print(key + ", " + value + "\n")

    title_init = "Histogram - Degree To Vtx - Initial Graph - Httpd_df"
    xlabel_init = "outDegree"
    ylabel_init = "vtxCnt"
    filename_init = "histogram-degree-to-vtx-init.pdf"

    if os.path.exists(filename_init):
        os.remove(filename_init)

    generate_histo(degree_to_vtx_dict_init, max_degree_init, title_init, xlabel_init, ylabel_init, 10, filename_init)

    title_final = "Histogram - Degree To Vtx - Final Graph - Httpd_df"
    xlabel_final = "outDegree"
    ylabel_final = "vtxCnt"
    filename_final = "histogram-degree-to-vtx-final.pdf"

    if os.path.exists(filename_final):
        os.remove(filename_final)

    generate_histo(degree_to_vtx_dict_final, max_degree_final,  title_final, xlabel_final, ylabel_final, 10, filename_final)


def main():
    # Replace 'yourfile.txt' with the path to the file you want to read
    file_path = 'httpd_df.txt'
    
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
            process_out_degree_to_vtx(lines)

    except FileNotFoundError:
        print("File not found.")
    except IOError:
        print("Error while reading the file.")

if __name__ == "__main__":
    main()