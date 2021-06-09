import numpy as np
import matplotlib.pyplot as plt


def stacked_bar(data, series_labels=None, category_labels=None,
                show_values=False, value_format="{}", y_label=None,
                grid=True, reverse=False, y_limit=None, size_plot=None, use_dataframe=False, throw_zeros=False,dict_colors={}):
    """Plots a stacked bar chart with the data and labels provided.

    Keyword arguments:
    data            -- 2-dimensional numpy array or nested list containing data for each series in rows
    series_labels   -- list of series labels (these appear in the legend)
    category_labels -- list of category labels (these appear on the x-axis)
    show_values     -- If True then numeric value labels will be shown on each bar
    value_format    -- Format string for numeric value labels (default is "{}")
    y_label         -- Label for y-axis (str)
    grid            -- If True display grid
    reverse         -- If True reverse the order that the series are displayed (left-to-right or right-to-left)
    y_limit         -- containes a int\float that will be the highest y value shown in the graph and y axis
    size_plot       -- contains an array of [ width , hight] we want the plot square area size will be
    use_dataframe   -- Bool, if true, data is treated as pandas df with series labels and category labels as rows and colums respectivly
    throw_zeros     -- Only applicable if use_dataframe is True, throws rows with all zeros in them
    """
    if throw_zeros and not use_dataframe:
        # TODO make throw zeros work without df too
        raise ValueError("throw_zeros only works if use_dataframe is chosen")

    # if throw zeros, remove rows with all zeros
    if throw_zeros:
        data = data[(data.T != 0).any()]
    # if data frame extract info from dataframe
    if use_dataframe:
        # remove no_change filter if needed:
        if 'no_change' in data.index:
            data = data.drop(['no_change'])
        series_labels = data.index
        category_labels = data.columns
        data = data.values
    ny = len(data[0])

    ind2 = range(ny)

    axes = []
    cum_size = np.zeros(ny)

    data = np.array(data)

    if reverse:
        data = np.flip(data, axis=1)
        category_labels = reversed(category_labels)

    if size_plot:
        fig = plt.figure(figsize=size_plot)
    plt.rcParams['font.size'] = '20'

    suit_colors_dict = {}
    for index, column in enumerate(series_labels):
        suit_colors_dict[index] = dict_colors[column]
    #print(data)
    sum_column = np.sum(data, axis=0)
    #print("old_data",data)
    #print("sum_column", sum_column)
    data = data.astype(float)
    for row_index in range(len(data)):
        for column_index in range(len(data[row_index])):
            if data[row_index][column_index] != 0.0:
                #print("before", "data[row_index][column_index]",data[row_index][column_index],"sum_column[column_index]*100", sum_column[column_index]*100)
                data[row_index][column_index] = format(data[row_index][column_index]/sum_column[column_index]*100, '.2f')
                #print("after:","\n","data[row_index][column_index]",data[row_index][column_index])
    #print("new data", data)
    #print("category_labels",category_labels )
    #print("series_labels",series_labels)
    # set the text in the same color as the bar
    for i, row_data in enumerate(data):
        axes.append(plt.bar(ind2, row_data, bottom=cum_size,
                            label=series_labels[i]))
        for row in range(len(row_data)):
            axes[i][row].set_color(suit_colors_dict[i])
        cum_size += row_data

    if not category_labels is None:
        plt.xticks(ind2, category_labels, rotation=20, fontsize=30)

    if y_label != None:
        plt.ylabel(y_label, fontsize=30)

    plt.legend()

    if grid:
        plt.grid()

    if y_limit != None:
        plt.ylim(0, y_limit)

    if show_values:
        max_tmp = []
        for axis in axes:
            max_tmp.append(max([bar.get_height() for bar in axis]))
        max_height_data = max(max_tmp)
        proportion_to_high = 0.08*max_height_data
        need_arrow = 0.08*max_height_data
        start_extra_heights = [axes[-1][i].get_y() + axes[-1][i].get_height() for i in range(len(axes[-1]))]
        jumps = [proportion_to_high for i in range(len(axes[0]))]
        for index,axis in enumerate(axes):
            for counter, bar in enumerate(axis):
                max_height = start_extra_heights[counter]
                w, h = bar.get_width(), bar.get_height()
                if 0.0 < h < need_arrow:
                    plt.annotate(value_format.format(h)+'%', xy=(bar.get_x(), bar.get_y()),
                                 xytext=(bar.get_x() + 0.2, max_height + jumps[counter]), color=suit_colors_dict[index],
                                 arrowprops=dict(arrowstyle="->"))
                    jumps[counter] += proportion_to_high * 1.2
                elif h > 0.0:
                    plt.text(bar.get_x() + w / 2, bar.get_y() + h / 2, value_format.format(h)+'%', ha="center",
                             va="center")
        # adding the number of total lines of the original pileups
        for index, bar in enumerate(axes[-1]):
            max_height = start_extra_heights[index]
            if max_height == 0.0:
                max_height = 1.3
            plt.annotate(value_format.format(sum_column[index]), xy=(bar.get_x(), bar.get_y()+bar.get_height()),
                         xytext=(bar.get_x(), max_height + jumps[index]),
                         arrowprops=dict(arrowstyle='fancy'))

    return plt, axes
