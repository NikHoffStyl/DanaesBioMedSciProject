import math
import os
import pandas as pd
import matplotlib.pyplot as plt


def convertRatioToFloat(x):
    a,b = x.split("/")
    c = float(int(a)/int(b))
    return c

def getLoge(x):
    lnx = math.log(x)
    return lnx

def getDF(fileName, limPvalue=None, verbose=False):
    df = pd.read_csv("data/"+fileName, sep='\t', header=0)
    df['Overlap'] = df['Overlap'].apply(convertRatioToFloat)
    df['Log P-value'] = df['P-value'].apply(getLoge)
    df['Log Adjusted P-value'] = df['Adjusted P-value'].apply(getLoge)
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

    df.drop(["Old P-value", "Old Adjusted P-value"], axis=1, inplace=True)

    return df

def getPathCountPerGene(dataL, genelist=None, limitPvalue=None):
    pathCountList = []
    diction = {}

    for fIndex, fName in enumerate(dataL):
        df = getDF(fName, limitPvalue)
        for index in df.index:
            geneSubList = df.loc[index, 'Genes']
            geneSubList = geneSubList.split(";")
            for ig, gene in enumerate(geneSubList):
                if fName+"__"+gene not in diction.keys():
                    diction.update({fName+"__"+gene : 1})
                else:
                    nPaths = diction[fName+"__"+gene] + 1
                    diction.update({fName+"__"+gene: nPaths})
        del df

    gList = [x.split("__")[1] for x in diction.keys()]
    gList = list(dict.fromkeys(gList))
    if genelist is not None: print("Gene List Comparison: ", genelist == gList)
    fList = [x.split("__")[0] for x in diction.keys()]
    fList = list(dict.fromkeys(fList))
    if dataL is not None: print("File List Comparison: ", dataL == fList)


    return diction, gList, pathCountList, fList

def mkdirIfNotExists(outDirectory_, opts=''):
    if not os.path.exists(outDirectory_):
        # print(os.path)
        print(os.getcwd())
        os.system('mkdir ' + opts + ' ' + outDirectory_)
    # print(os.path)
    print(os.getcwd())

def main(indexOfFile=None, numberOfFiles=None, compact=False, combCellBioMol=False, limitPvalue=None):
    dataList = [x for x in os.walk("data").__next__()[2] if ".txt" in x and "KEGG" not in x]
    colorList = ['maroon', 'springgreen', 'royalblue', 'violet', 'gold', 'silver',
                 'red', 'green','deepskyblue', 'mediumvioletred', 'darkorange', 'sienna']
    if numberOfFiles is not None:
        if indexOfFile is not None: return -1
        else:
            dataList = dataList[:numberOfFiles]
            colorList = colorList[:numberOfFiles]
            imageName = "total_"
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
            imageName = "imagesPerFile"
        if limitPvalue is not None: imageName += "lim" + str(limitPvalue).replace("0.0", "")
        mkdirIfNotExists(imageName)
        imageName += "/"+str(indexOfFile)+"_"
    else:
        imageName = "total__"
        if limitPvalue is not None: imageName += "lim" + str(limitPvalue).replace("0.0", "") + "_"

    print(dataList)
    print("Number of Files Conbsidered: %d" % len(dataList))
    diction1, g_list, p_list, f_list = getPathCountPerGene(dataList, limitPvalue=limitPvalue)

    stepDF = pd.DataFrame(columns=dataList+["Total Count"],index=g_list)
    for gene in g_list:
        seriesPerGene = {}
        totalCount = 0
        for fileN in f_list:
            if fileN+"__"+gene in diction1:
                seriesPerGene.update({fileN:diction1[fileN+"__"+gene]})
                totalCount += diction1[fileN+"__"+gene]
            else:
                seriesPerGene.update({fileN: 0})
        seriesPerGene.update({"Total Count": totalCount})
        stepDF.loc[gene] = pd.Series(seriesPerGene)
    stepDF.sort_values(["Total Count"], inplace=True)
    # print(stepDF.info())
    # print(stepDF.head().to_string(justify="left"), end="\n\n")
    # print(stepDF.tail().to_string(justify="left"), end="\n\n")
    stepDF.drop(["Total Count"], axis=1, inplace=True)
    numGenes = len(g_list)

    if compact:
        for ii, gn in enumerate(stepDF.index):
            if ii <(numGenes - 40): stepDF.drop(gn, inplace=True)
        stepDF.plot.barh(figsize=(10,10), stacked=True, color=colorList)
        imageName += "compactGeneList.png"
    else:
        if indexOfFile is None: stepDF.plot.barh(figsize=(10,18), stacked=True, color=colorList)
        if numberOfFiles is None: stepDF.plot.barh(figsize=(10,200), stacked=True, color=colorList)
        imageName += "fGeneList.png"
    plt.margins(y=0)
    plt.savefig(imageName, bbox_inches='tight')
    plt.show()

    plt.close('all')

    return 0

def drawHist(datList, dfColumn, xLabel, dirName, colour, limitPvalue):
    for fIndex, fName in enumerate(datList):
        datFrame = getDF(fName, limPvalue=limitPvalue)
        datFrame[dfColumn].plot(kind='hist', bins=25, alpha=0.75, color=[colour[fIndex]])
        plt.xlabel(xLabel)
        plt.ylabel('Number of Pathways')
        plt.margins(y=0.1)
        pngName = dirName + dfColumn.replace(" ","")
        plt.savefig(pngName, bbox_inches='tight')
        plt.show()
    plt.close('all')
    del datFrame

    return True

def drawScatter(datList, dfColumn, dfColumn2, dirName, colour, limitPvalue):
    for fIndex, fName in enumerate(datList):
        datFrame = getDF(fName, limPvalue=limitPvalue)
        datFrame.plot(kind='scatter', x=dfColumn, y=dfColumn2, alpha=0.75, color=[colour[fIndex]])
        yLabel = dfColumn2.replace("Log P-value", "ln(p-value)")
        plt.ylabel(yLabel)
        plt.margins(y=0.1)
        pngName = dirName + dfColumn.replace(" ","") + "Vs" + dfColumn2.replace(" ","")
        plt.savefig(pngName, bbox_inches='tight')
        plt.show()
    plt.close('all')
    del datFrame

    return True

def plotAll(indexOfFile=None, numberOfFiles=None, combCellBioMol=False, limitPval=None):
    dataList = [x for x in os.walk("data").__next__()[2] if ".txt" in x and "KEGG" not in x]
    colorList = ['maroon', 'springgreen', 'royalblue', 'violet', 'gold', 'silver',
                 'red', 'green','deepskyblue', 'mediumvioletred', 'darkorange', 'sienna']
    if numberOfFiles is not None:
        if indexOfFile is not None: return -1
        else:
            dataList = dataList[:numberOfFiles]
            colorList = colorList[:numberOfFiles]
            imageName = "total_"
            if limitPval is not None: imageName += "lim" + str(limitPval).replace("0.0","") + "_"
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
            imageName = "imagesPerFile"
        if limitPval is not None: imageName += "lim" + str(limitPval).replace("0.0","")
        mkdirIfNotExists(imageName)
        imageName += "/"+str(indexOfFile)+"_"
    else:
        imageName = "total__"
        if limitPval is not None: imageName += "lim" + str(limitPval).replace("0.0","") + "_"

    print(dataList)
    print("Number of Files Conbsidered: %d" % len(dataList))


    mkdirIfNotExists(imageName + "histograms", opts="-p")

    drawHist(dataList, 'Overlap', 'Overlap', imageName+'histograms/', colorList, limitPval)
    drawHist(dataList, 'Log P-value', 'ln(p-value)', imageName+'histograms/', colorList, limitPval)
    drawHist(dataList, 'Log Adjusted P-value', 'ln(Adjusted p-value)', imageName+'histograms/', colorList, limitPval)
    drawHist(dataList, 'Odds Ratio', 'Odds Ratio', imageName+'histograms/', colorList, limitPval)
    drawHist(dataList, 'Combined Score', 'Combined Score', imageName+'histograms/', colorList, limitPval)
    drawHist(dataList, 'P-value', 'p-value', imageName+'histograms/', colorList, limitPval)

    mkdirIfNotExists(imageName + "scatter", opts="-p")
    drawScatter(dataList, 'Combined Score', 'Log P-value', imageName + "scatter/", colorList, limitPval)
    drawScatter(dataList, 'Odds Ratio', 'Log P-value', imageName + "scatter/", colorList, limitPval)
    drawScatter(dataList, 'Overlap', 'Log P-value', imageName + "scatter/", colorList, limitPval)
    drawScatter(dataList, 'Log Adjusted P-value', 'Log P-value', imageName + "scatter/", colorList, limitPval)
    drawScatter(dataList, 'Adjusted P-value', 'P-value', imageName + "scatter/", colorList, limitPval)

    return True


if __name__ == '__main__':
    for i in range(12):
        main(indexOfFile=i, limitPvalue=0.05, compact=True)
        main(indexOfFile=i, limitPvalue=0.05, compact=True, combCellBioMol=True)

        main(indexOfFile=i, limitPvalue=0.05)
        main(indexOfFile=i, limitPvalue=0.05, combCellBioMol=True)

        main(indexOfFile=i, compact=True)
        main(indexOfFile=i, compact=True, combCellBioMol=True)

        main(indexOfFile=i)
        main(indexOfFile=i, combCellBioMol=True)

    main(limitPvalue=0.05,compact=True)
    main(limitPvalue=0.05)
    main(compact=True)
    main()

    for i in range(12):
        plotAll(indexOfFile=i)
        plotAll(indexOfFile=i, limitPval=0.05)
        plotAll(indexOfFile=i, combCellBioMol=True)
        plotAll(indexOfFile=i, combCellBioMol=True, limitPval=0.05)

    print("Finished Running!")
