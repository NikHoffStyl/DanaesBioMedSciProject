# import math
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import tools


def getDF(fileName, limPvalue=None, droplist=None, verbose=False):
    df = pd.read_csv("data/"+fileName, sep='\t', header=0)
    df['Overlap'] = df['Overlap'].apply(tools.convertRatioToFloat)
    df['Log P-value'] = df['P-value'].apply(tools.getLoge)
    df['Log Adjusted P-value'] = df['Adjusted P-value'].apply(tools.getLoge)
    df['Gene Count'] = df['Genes'].apply(tools.countGenes)
    df['Pathway'] = df['Term'].apply(tools.rmGOkey)
    if verbose: print(df.info())
    keywordList = [" bone", " eye", " extrem", " nose", " beak", " tail", " claw", " ear"]

    for index in df.index:
        if any([df.loc[index, 'Term'].find(keyword) != -1 for keyword in keywordList]):
            if verbose: print(index, df.loc[index, 'Term'])
            df.drop(index, inplace=True)
        else:
            if limPvalue is not None:
                if df.loc[index, 'P-value'] >= limPvalue:
                    df.drop(index, inplace=True)
    if droplist is None:
        df.drop(["Old P-value", "Old Adjusted P-value", "Log Adjusted P-value","P-value", "Adjusted P-value",
                "Overlap", "Odds Ratio", "Combined Score", "Term"], axis=1, inplace=True)
    else:
        df.drop(droplist, axis=1, inplace=True)

    return df


def main(indexOfFile=None, numberOfFiles=None, compact=False, combCellBioMol=False, limitPvalue=None):
    dataList = [x for x in os.walk("data").__next__()[2] if ".txt" in x and "KEGG" not in x]
    colorList = ['maroon', 'springgreen', 'royalblue', 'violet', 'gold', 'silver',
                 'red', 'green','deepskyblue', 'mediumvioletred', 'darkorange', 'sienna']
    if numberOfFiles is not None:
        if indexOfFile is not None: return -1
        else:
            dataList = dataList[:numberOfFiles]
            colorList = colorList[:numberOfFiles]
            imageName = "imagesMaxPval2GeneCount"
            tools.mkdirIfNotExists(imageName)
            imageName += "/total_"
            if limitPvalue is not None: imageName += "lim" + str(limitPvalue).replace("0.0", "") + "_"
    elif indexOfFile is not None:
        if combCellBioMol:
            if "molecular" not in dataList[indexOfFile]: return -1
            else:
                dataList_bkup = dataList.copy()
                dataList = [dataList[indexOfFile], dataList[indexOfFile].replace("molecular", "cell"),
                            dataList[indexOfFile].replace("molecular", "bio")]
                colorList = [colorList[dataList_bkup.index(x)] for x in dataList]
                imageName = "imagesBioCellMol"
        else:
            dataList = [dataList[indexOfFile]]
            colorList = [colorList[indexOfFile]]
            imageName = "imagesMaxPval2GeneCount/perFile"
        if limitPvalue is not None: imageName += "lim" + str(limitPvalue).replace("0.0", "")
        tools.mkdirIfNotExists(imageName, opts="-p")
        imageName += "/"+str(indexOfFile)+"_"
    else:
        imageName = "total__"
        if limitPvalue is not None: imageName += "lim" + str(limitPvalue).replace("0.0", "") + "_"

    print(dataList)
    print("Number of Files Conbsidered: %d" % len(dataList))

    # # fig, axes = plt.subplots(nrows=1, ncols=2)
    # fig, ax = plt.subplots()
    for fName in dataList:
        df = getDF(fName, limitPvalue)
        df.sort_values(["Log P-value"], inplace=True)

        nRows = len(df.index)
        if compact:
            for ii, gn in enumerate(df.index):
                if ii >10: df.drop(gn, inplace=True)

            DFumpy_0 = df["Pathway"].to_numpy()
            DFumpy_1 = df["Gene Count"].to_numpy()
            DFumpy_2 = df["Log P-value"].to_numpy()

            ax1 = plt.subplot(132)
            plt.barh(DFumpy_0, DFumpy_2, color=colorList[0])
            ax1.set(xlabel="P-value", ylabel="Pathways", title=fName.replace("_TABLE.txt", ""))

            ax2 = plt.subplot(133, sharey=ax1)
            plt.barh(DFumpy_0, DFumpy_1, color=colorList[0], alpha=0.5)
            plt.tick_params('y', labelleft=False)
            ax2.set(xlabel="Gene Count")

            imageName += "compact.png"
        # else:
        #     if indexOfFile is None: df.plot.barh(figsize=(10,18), stacked=True, color=colorList)
        #     if numberOfFiles is None: df.plot.barh(figsize=(10,200), stacked=True, color=colorList)
        #     imageName += "full.png"
        plt.margins(y=0)
        plt.savefig(imageName, bbox_inches='tight')
        plt.show()

        plt.close('all')

    return 0




if __name__ == '__main__':
    for i in range(12):
        main(indexOfFile=i, limitPvalue=0.05, compact=True)

    # main(limitPvalue=0.05,compact=True)
    # main(limitPvalue=0.05)
    # main(compact=True)
    # main()


    print("Finished Running!")
