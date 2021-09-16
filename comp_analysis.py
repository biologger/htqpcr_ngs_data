#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import seaborn as sns
import math
import json
import itertools
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.colors import ListedColormap
from scipy.stats import linregress
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib.ticker import FuncFormatter


class MathFunctions:

    @staticmethod
    def calc_rel_abund(df, selection=None):
        if isinstance(selection, list):
            return 100/df.loc[selection].sum()*df.loc[selection]
        if selection is None:
            return 100/df.sum()*df
        print("Selection has to be a list")

    @staticmethod
    def center_elts(x):
        return np.exp(np.log(x) - np.mean(np.log(x)))

    @staticmethod
    def log10_tf(df, col):
        return df[col].apply(np.log10).replace([-np.inf, np.inf], 0)

    @staticmethod
    def bray_curtis_dissimilarity(A, B):
        u, v = np.array(A), np.array(B)
        return abs(u-v).sum()/abs(u+v).sum()


class HelperFunctions:


    def newtax_dict(self):
        # dictionary for new taxonomic classification of Lactobaciilus group
        with open('rawdata/new_species_names_dict.json') as f:
            for line in f:
                ntax_dict = json.loads(line)
        return ntax_dict

    def find_16S_copy_numbers(self, inputfile, specieslist, average_16S_CN):
        ntax_dict = self.newtax_dict()
        rev_dict = {v: k for k, v in ntax_dict.items()}
        copynum = pd.read_csv(inputfile, sep="\t")
        copy_sel = copynum.set_index('name')
        copy_data = copy_sel[copy_sel.index.str.contains("|".join(specieslist))]
        missing = set(specieslist) - set(copy_data.index)
        old_names = [rev_dict[o] for o in missing if o in rev_dict.keys()]
        old_data = copy_sel[copy_sel.index.str.contains("|".join(old_names))]
        missing = missing - set([ntax_dict[n] for n in old_data.index])
        h_taxa = list(set([o.split(" ")[0] for o in missing])) + ["Lactobacillus"]
        taxa_data = copy_sel[
            (copy_sel.index.str.contains("|".join(h_taxa))) &
            (copy_sel["rank"] != "species")]
        missing = set(h_taxa) - set(taxa_data.index)
        vals = [[np.nan]*len(copy_sel.columns) for m in missing]
        missing_data = pd.DataFrame(vals, index=missing, columns=copy_sel.columns)
        cn_df = pd.concat([
            copy_data.sort_index(), old_data.sort_index(),
            taxa_data.sort_index(), missing_data.sort_index()], axis=0)
        tbl_spec_names = [ntax_dict[s] if s in ntax_dict.keys() else s for s in list(cn_df.index)]
        cn_df.index = tbl_spec_names
        copy_tbl = cn_df[["rank", "childcount", "min", "max", "median", "mean", "stddev"]]
        copy_tbl.columns = ['Rank', 'N', 'Min.', 'Max.', 'Median', 'Avg.', 'SD']
        copy_tbl.to_csv("rawdata/16S_copy_number_table.csv")

        # create a 16S copy number dictionary
        copy_dict = dict(copy_tbl["Avg."])
        # Include L. plantarum group
        lplantgroup = copy_tbl[copy_tbl.index.str.contains("plantarum|pentosus")]
        lplantarumgroup_m = lplantgroup["Avg."].mean().round()
        copy_dict.update({"L. plantarum group": lplantarumgroup_m})
        # add info from higer taxonomic rank
        missing = set(specieslist) - set(copy_dict.keys())
        for l in missing:
            data = copy_tbl.loc[l.split(" ")[0], "Avg."]
            copy_dict.update({l: data})
            if np.isnan(data):
                print("Missing data for {}, set average CN ({})".format(l, average_16S_CN))
                copy_dict.update({l: average_16S_CN})

        with open("rawdata/copydict.json", "w") as f:
            f.write(json.dumps(copy_dict))

        return copy_tbl, copy_dict


    def get_assay_dict(self):
        map_file = os.path.join("HTqPCR_dataparser", "labelfiles", "assay_species.csv")
        adf = pd.read_csv(map_file, header=None)
        assay_dict = dict(zip(adf[0].values, adf[1].values))
        return assay_dict


    def new_assay_species_labels(self, specieslist):
        ntax_dict = self.newtax_dict()
        assay_dict = self.get_assay_dict()
        newlabel = []
        for label in specieslist:
            newl = assay_dict[label]
            if newl in ntax_dict:
                newl = ntax_dict[newl]
            newlabel.append(newl)
        return newlabel


    def new_species_labels(self, specieslist):
        ntax_dict = self.newtax_dict()
        newlabel = []
        for label in specieslist:
            newl = " ".join(label.split("_"))
            if newl in ntax_dict:
                newl = ntax_dict[newl]
                if 'subsp.' in newl:
                    newl = newl.split(' subsp.')[0]
            newlabel.append(newl)
        return newlabel


    @staticmethod
    def sort_sum(df, ascending=False):
        df["sort"] = df.sum(axis=1)
        df.sort_values("sort", ascending=ascending, inplace=True)
        df.drop("sort", axis=1, inplace=True)
        return df


    @staticmethod
    def get_stat_df(data_df):
        data_stat = pd.DataFrame()
        data_df = data_df.replace(0, np.nan)
        data_stat["Mean"] = data_df.T.mean()
        data_stat["Std"] = data_df.T.std()
        data_stat["Min"] = data_df.T.min()
        data_stat["Max"] = data_df.T.max()
        data_stat["Median"] = data_df.T.median()
        data_stat["Count"] = (data_df.T > 0).sum()
        data_stat.sort_values(by=["Count"], ascending=False, inplace=True)
        return data_stat


    @staticmethod
    def create_summarydf(qpcr_count, ngs_count, qpcr_reldna, ngs_reldna):
        cheese_data = pd.DataFrame({
            "qPCR_count": qpcr_count.unstack(), "NGS_count": ngs_count.unstack(),
            "qPCR_rel": qpcr_reldna.unstack(), "NGS_rel": ngs_reldna.unstack()})
        cheese_data = cheese_data.replace(np.nan, 0)
        shared_species = list(qpcr_count.index.intersection(ngs_count.index))

        cheese_dat = cheese_data.reset_index()
        col_names = ["Sample", "Species"] + list(cheese_data.columns)
        cheese_dat.columns = col_names
        cheese_dat["shared_positive"] = (cheese_dat["NGS_count"] != 0) & (cheese_dat["qPCR_count"] != 0)

        cheese_dat["qPCR_only"] = (cheese_dat["qPCR_count"] != 0) & (cheese_dat["NGS_count"] == 0)
        cheese_dat["NGS_only"] = (
            (cheese_dat.reset_index()["Species"].str.contains("|".join(shared_species)))
            & (cheese_dat["NGS_count"] != 0) & (cheese_dat["qPCR_count"] == 0))

        cheese_dat["NGS_exclusive"] = (
            (~cheese_dat.reset_index()["Species"].str.contains("|".join(shared_species)))
            & (cheese_dat["NGS_count"] > 0))

        cheese_dat["not_detected"] = (cheese_dat["NGS_count"] == 0) & (cheese_dat["qPCR_count"] == 0)

        cheese = cheese_dat.melt(id_vars=col_names,
                var_name="Category",
                value_name="Cat_value")

        cheese_df = cheese[cheese["Cat_value"] == True]
        cheese_df = cheese_df.reset_index(drop=True)
        cheese_df = cheese_df.drop("Cat_value", axis=1)
        return cheese_df, shared_species


    def reads_info(df, verbose=True):
        if verbose:
            # NGS analysis stats
            seqdepth = int(df.iloc[:, 9::].sum().mean().round())
            seqd_min = df.iloc[:, 9::].sum().min().round()
            seqd_max = df.iloc[:, 9::].sum().max().round()
            totalreads = df.iloc[:, 9::].sum().sum()
            print("Average sequencing depth: {} reads\nRange: {}-{} reads".format(
                seqdepth, seqd_min, seqd_max))
            print("Total assigned reads: {}\nNumber of ASVs: {}".format(
                totalreads, len(df.index)))

    @staticmethod
    def tax_info(df, verbose=True, printlimit=100):
        # Taxonomic information
        if verbose:
            for c in list(df.iloc[:, 2:9].columns.values):
                diff = list(set(df[c].values))
                identified = []
                no_spec_id = []
                for d in diff:
                    if c in d:
                        no_spec_id.append(d)
                    else:
                        identified.append(d)
                print("\n{} different groups found in {}\n".format(len(diff), c))
                if len(diff) < printlimit:
                    print(diff)
                if len(no_spec_id) > 0:
                    print("\n{} identified, {} not identified: {}".format(
                        len(identified), len(no_spec_id), ", ".join(no_spec_id)))


    def asv_to_speciestable(self, df):
        # Combine ASVs to species
        species = list(set(df["Species"]))
        spec_dict = {}
        for spec in species:
            q = df[df["Species"] == spec]
            q.iloc[:, 9:]
            spec_dict.update({spec: q.iloc[:, 9:].sum()})

        ngs_count = pd.DataFrame.from_dict(spec_dict).T
        ngs_count.index = self.new_species_labels(ngs_count.index)
        return ngs_count


    def method_comparison_df(methods):
        n_smpl = len(methods["R"].columns)
        df = pd.DataFrame()
        for k in methods.keys():
            col = []
            for i in range(n_smpl):
                col.append(k + str(i+1))
            methods[k].columns = col
            df = df.append(methods[k].T, sort=True)

        sorted_cols = [itm + str(i+1) for i in range(n_smpl) for itm in methods.keys()]
        df_r = df.T.replace(np.nan, 0)
        df_r = df_r[sorted_cols]
        return df_r


    def dissimilarity_df(df, methods=["R", "A", "B"], dist_calc='braycurtis'):
        data = []
        n_met = len(methods)
        for j in range(1, n_met):
            k=0
            for i in range(0, len(df.columns), n_met):
                R, X = df.iloc[:, i], df.iloc[:, i+j]
                if dist_calc == 'braycurtis':
                    diss = MathFunctions().bray_curtis_dissimilarity(R, X)
                data.append(["S{:02}".format(k+1), methods[j], diss])
                k += 1
        diss_df = pd.DataFrame(data, columns=["Sample", "Method", "Dissimilarity"])
        return diss_df

    @staticmethod
    # Labels with abbreviatet species names
    def abbrev_species_labels(specieslabels):
        return [x.split(" ")[0][0] +". " + " ".join(x.split(" ")[1:]) for x in specieslabels]


    @staticmethod
    # get labels with number of observations (Figure 3)
    def get_better_labels(df, labels, col_name="Category"):
        newlabels = []
        for label in labels:
            if label == col_name:
                new_label = label
            else:
                number = len(df[df[col_name] == label])
                new_label = "{} ({})".format(" ".join(label.split("_")), number)
            newlabels.append(new_label)
        return newlabels


    def get_xrange(min_val, max_val):
        logscale_min = math.floor(np.log10(min_val))
        logscale_max = math.ceil(np.log10(max_val))
        xlim = (10**logscale_min, 10**logscale_max)
        ticks = [10**(i) for i in np.arange(logscale_min, logscale_max+1, 1, dtype=float)]
        return xlim, ticks


class PlotFunctions:

    def draw_qpcr_heatmap(qpcrdata, annotation, ax, cax):
            cmap = ListedColormap(sns.color_palette("Blues", 5))
            cmap.set_under("lightgrey")
            cbar_label = "log copies/\u03BCl"

            sns.heatmap(
                qpcrdata, vmin=3, vmax=8, cmap=cmap, annot=annotation,
                annot_kws={"size": 8, "color": "black"}, fmt='',
                cbar_kws={
                        'orientation': 'vertical', "label": cbar_label,
                        'extend': "min"},
                ax=ax, cbar_ax=cax, linewidths=1, linecolor="black")

            newlabels = []
            for label in ax.get_ymajorticklabels():
                if "subsp." in label.get_text():
                    italic = label.get_text().split("subsp.")
                    newtext = (
                        '$\it{' + "\ ".join(italic[0].split(" ")) + '}$'
                        + "subsp. " + '$\it{' + italic[1].strip() + '}$')
                else:
                    newtext = '$\it{' + "\ ".join(label.get_text().split(" ")) + '}$'

                label.set_text(newtext)
                newlabels.append(label)

            ax.set_yticklabels(newlabels, rotation=0, fontsize=12)
            ax.set_xticklabels(ax.get_xmajorticklabels(), rotation=0, fontsize=12)
            ax.tick_params(axis='both', length=0.0, width=0.0)
            ax.yaxis.tick_right()
            ax.yaxis.set_label_position('right')

            return ax


    def draw_ngs_barplot(df_high, df_low, ax, th):
        l_palette = sns.color_palette("Paired") + [sns.color_palette("Greys",5)[2]]
        # Upper plot > th% rel. abundance
        q = df_high.T.plot.bar(
                stacked=True, ax=ax[0], legend=False, sharex=False,
                title="≥ " + str(th) + "% average abundance")
        # Create first legend
        handles, labels = q.get_legend_handles_labels()
        first_legend = plt.legend(
                            handles, labels, loc='upper right',
                            bbox_to_anchor=(1.40, 2.146),
                            frameon=False, prop={"style": "italic"})
        plt.gca().add_artist(first_legend)
        # Lower plot < th% rel abundance
        q2 = df_low.T.plot.bar(
                stacked=True, ax=ax[1], legend=False, color=l_palette,
                sharex=False, title="< " + str(th) + "% average abundance")
        # Legend
        handles, labels = q2.get_legend_handles_labels()
        leg = plt.legend(
            handles, labels, loc='upper right', frameon=False,
            prop={"style": "italic"}, bbox_to_anchor=(1.408, 1.023))
        # y-axis labels
        ax[0].set_ylabel("relative abundance [%]")
        ax[1].set_ylabel("relative abundance [%]")
        # remove ticks
        ax[0].tick_params(axis='x', length=0.0, width=0.0, which="both", rotation=0)
        ax[1].tick_params(axis='x', length=0.0, width=0.0, which="both", rotation=0)

        for txt in leg.get_texts():
            if "Other species" in txt.get_text():
                txt.set_style("normal")

        return ax


    def count_data_plot(df, ax, color_pal):
        # linear reg calculations
        df["NGS_log"] = MathFunctions().log10_tf(df, "NGS_count")
        df["qPCR_log"] = MathFunctions().log10_tf(df, "qPCR_count")
        shared = df.query('Category == "shared_positive"')
        slope, intercept, r_value, p_value, std_err = linregress(
            shared["NGS_log"], shared["qPCR_log"])
        colors = color_pal[0: len(set(df["Category"]))]
        sns.scatterplot(
            x="NGS_log", y="qPCR_log", data=df,
            hue="Category", palette=colors, ax=ax, s=60, alpha=0.8)
        if p_value <= 0.05:
            ax.annotate(
                "R² = {:.3f}".format(r_value ** 2),
                xy=(0.025, 0.9), fontsize=10, xycoords='axes fraction')
        ax.set_ylabel("HT-qPCR\nlog(copies/\u03BCl)", fontsize=10)
        ax.set_xlabel("log(reads)\nNGS", fontsize=10)
        ax.axhline(np.log10(800), color='red', linestyle='--', alpha=0.7)

        # legend
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(
                        handles,
                        HelperFunctions().get_better_labels(df, labels),
                        ncol=2, frameon=False,
                        loc=4, bbox_to_anchor = (1.01, 0.1),
                        columnspacing=0.4, labelspacing=0.5,
                        handletextpad=0.25, prop={"size": 10})
        return ax

    def rel_data_plot(df, ax, color_pal, mode="rel"):
        colors = color_pal[0: len(set(df["Category"]))]
        df = df.query('Category != "NGS_exclusive"')
        colors.pop(len(set(df["Category"]))-1)

        if mode != "rel":
            df["qPCR_log"] = MathFunctions().log10_tf(df, "qPCR_count")
            yd = "qPCR_log"
            ax.axhline(np.log10(800), color='red', linestyle='--', alpha=0.7)
        else:
            yd = "qPCR_rel"
            shared = df.query('Category == "shared_positive"')
            slope, intercept, r_value, p_value, std_err = linregress(
                shared["NGS_rel"], shared["qPCR_rel"])
        sns.scatterplot(
            x="NGS_rel", y=yd, data=df,
            hue="Category", palette=colors, ax=ax, s=60, alpha=0.9)
        if mode != "rel":
            ax.set_ylabel("HT-qPCR\nlog(copies/\u03BCl)", fontsize=10)
        else:
            ax.set_ylabel("HT-qPCR\nrelative abundance [%]", fontsize=10)
            if p_value <= 0.05:
                ax.annotate(
                    "R² = {:.3f}".format(r_value ** 2),
                    xy=(0.025, 0.9), fontsize=10, xycoords='axes fraction')
        ax.set_xlabel(
            "relative abundance [%]\nNGS",
            fontsize=10)

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(
            handles=handles,
            labels=HelperFunctions().get_better_labels(df, labels),
            loc=4, frameon=False, bbox_to_anchor = (1, 0.0),
            columnspacing=0.4, labelspacing=0.5, handletextpad=0.25,
            prop={"size": 10})

        return ax

    def dissimilarity_plot(df, ax, color_pal):
        df1 = df.query('Method == "A" or Method == "B"')
        df2 = df.query('Method == "C" or Method == "D"')

        clabels=["raw", "raw GCN", "corr.", "corr. GCN"]
        markers = ["o", "o", "^", "^"]
        cpal = color_pal[0:len(clabels)]
        leg_handle = []

        c1 = cpal[0:2]
        c2 = cpal[2::]
        sns.stripplot(
                x="Sample", y="Dissimilarity", hue="Method", data=df1,
                palette=c1, jitter=0., dodge=True, ax=ax, marker="o")
        sns.stripplot(
                x="Sample", y="Dissimilarity", hue="Method", data=df2,
                palette=c2, jitter=0., dodge=True, ax=ax, marker="^")

        for i, label in enumerate(clabels):
            handle = mlines.Line2D([], [], color=cpal[i], marker=markers[i],
                                   linestyle='None', label=label)
            leg_handle.append(handle)

        l1 = ax.legend(
                        handles=leg_handle, ncol=2, markerscale=1, framealpha=1,
                       columnspacing=0.4, labelspacing=0.5, handletextpad=0.25)

        ax.set_xlabel("")
        ax.set_ylabel(
            "Bray-Curtis dissmilarity",
            fontsize=10)
        ax.grid(True)
        ax.tick_params(axis='x', rotation=90)
        l1.get_frame().set_linewidth(0.0)

        return ax

    def linkage_plot(df, ax, shared_species, orient="left"):
        Z = linkage(df, method='average', metric="braycurtis", optimal_ordering=True)
        D = dendrogram(
            Z,
            orientation=orient,
            distance_sort=True,
            ax=ax, no_labels=True,
            color_threshold=0, above_threshold_color='black'
        )
        for spine in ax.spines:
            ax.spines[spine].set_visible(False)
        ax.set_xticklabels([])
        ax.tick_params(axis="x", width=0, length=0, pad=80)
        ax.set_ylabel("UPGMA linkage (Bray-Curtis dissimilarity)", fontsize=12)
        return ax, D["leaves"]

    def sort_barplots(df, leaves, shared_species):
        df = df.iloc[leaves,:]
        # Create a summary using shared and others
        mask = df.columns.isin(shared_species)
        spec_order = list(df[shared_species].mean(axis=0).sort_values(ascending=False).index)
        sumdf = df.loc[:, spec_order]
        sumdf["Other species"] = df.loc[:, ~mask].sum(axis=1)
        return sumdf

    def draw_linked_barplots(df, ax, legendax, refdict, cmap):
        # rel. abundance barplots
        q = df.plot.barh(
                stacked=True, ax=ax, legend=False, width=0.8,
                edgecolor='black', linewidth=0.8, color=cmap)
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()
        ax.set_xlim(0, 101)
        ax.tick_params(axis="y", width=0, length=0, pad=58)
        for spine in ax.spines:
            ax.spines[spine].set_visible(False)
        # labels
        labels = ["{} S{:02d}".format(refdict[l[0]], int(l[1:])) for l in list(df.index)]
        ax.set_yticklabels(labels, ha='right', fontsize=12)
        ax.tick_params(axis="x", labelsize=12)
        ax.set_xlabel("relative abundance [%]", fontsize=12)

        # Create legend
        handles, labels = q.get_legend_handles_labels()
        leg = legendax.legend(
                            handles, labels, loc=2, ncol=1,  bbox_to_anchor=(0.65, 1),
                            prop={"style": "italic", "size": 12}, frameon=False,
                            columnspacing=0.4, labelspacing=0.5, handletextpad=0.5,
                            borderpad=0.1)

        for txt in leg.get_texts():
            if "Other species" in txt.get_text():
                txt.set_style("normal")

        return ax

    @staticmethod
    def bias_plot(
            df1, df2=None, df3=None, df4=None, ax=None,
            color_pal="colorblind", legend=None, colors=None, markers=None):

        if ax is None:
            ax=plt.gca()
        datalist = []
        pos = np.arange(1, len(df1.index)+1,  dtype=float)
        sortdf = df1.sort_values(["estimate"], ascending=True)
        sorter = sortdf.index
        labels = HelperFunctions().abbrev_species_labels(
                                                    sortdf["taxon"].to_list())
        if df2 is None and df3 is None and df4 is None:
            df1.loc[sorter, "y_pos"] = pos
            datalist.append(df1)
        elif df3 is None and df4 is None:
            df1.loc[sorter, "y_pos"] = pos + 0.125
            df2.loc[sorter, "y_pos"] = pos - 0.125
            datalist.extend([df1, df2])
        elif df4 is None:
            df1.loc[sorter, "y_pos"]= pos + 0.125
            df2.loc[sorter, "y_pos"]  = pos
            df3.loc[sorter, "y_pos"]  = pos - 0.125
            datalist.extend([df1, df2, df3])
        else:
            df1.loc[sorter, "y_pos"]= pos + 0.25
            df2.loc[sorter, "y_pos"]  = pos + 0.125
            df3.loc[sorter, "y_pos"]  = pos - 0.125
            df4.loc[sorter, "y_pos"]  = pos - 0.25
            datalist.extend([df1, df2, df3, df4])

        if colors:
            colors=colors
        else:
            colors = sns.color_palette(color_pal, len(datalist))
        if markers:
            markers = markers
        else:
            markers = ["o", "o", "^", "^"]
        lim_x = []
        for i, data in enumerate(datalist):
            ax.plot(
                "estimate", "y_pos", marker=markers[i],
                color=colors[i], linewidth=0, data=data,
                markersize=6)
            ax.plot(
                [data["errorbar_min"], data["errorbar_max"]],
                [data["y_pos"], data["y_pos"]], color=colors[i], linewidth=2)

            lim_x.extend([data["errorbar_min"].min(), data["errorbar_max"].max()])

        # handle ticks
        ylim = ax.get_ylim()
        ax.plot([1, 1], list(ax.get_ylim()), color="black", linewidth=1)
        ax.set_yticks(pos)
        ax.set_yticklabels(labels, style="italic", fontsize=10)
        ylim = (ylim[0] + 0.5, ylim[1] - 0.5)
        ax.set_ylim(ylim)
        ax.set_xscale("log")
        min_x = min(lim_x) - 1/min(lim_x)
        max_x = max(lim_x) + 1/max(lim_x)
        ax.set_xlim(min_x, max_x)
        formatter = FuncFormatter(lambda x, _: '{:.16g}'.format(x))
        ax.xaxis.set_major_formatter(formatter)
        ax.grid(axis="x", which="major")

        ax.set_xlabel("Bias estimate", fontsize=10)

        return ax


def main():
    print("For usage in jupyter-notebooks or python scripts\nimport analysis_core")

if __name__ == "__main__":
    main()
