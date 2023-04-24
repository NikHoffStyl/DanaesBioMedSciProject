import os
import pandas as pd


def mkdirIfNotExists(outDirectory_, opts=''):
    if not os.path.exists(outDirectory_):
        os.system('mkdir ' + opts + ' ' + outDirectory_)


def main():
    df = pd.read_csv("data/raw/APPENDIX.csv", header=0)

    col_0  = df["BMPIinUP"].to_list()
    col_1 = df["BMPinDOWN"].to_list()
    col_2 = df["FGFinUP"].to_list()
    col_3 = df["FGFinUP_1"].to_list()
    col_4 = df["FGFinDOWN"].to_list()
    col_5 = df["FGFinDOWN_1"].to_list()
    col_6 = df["FGFinDOWN_2"].to_list()

    colSum = col_0 + col_1 + col_2 + col_3 + col_4 + col_5 + col_6
    colALl = list(dict.fromkeys(colSum))
    print("collAll:", colALl)
    for gene in colALl:
        if type(gene) != str :
            colALl.remove(gene)
            continue
    colALl.sort()
    print(colALl)
    colDict = dict.fromkeys(colSum)

    count = []
    for gene in colALl[:]:
        measure= colSum.count(gene)
        if measure <= 1:
            colALl.remove(gene)
            continue
        colDict[gene] = measure
        count.append(measure)

    newDF = pd.DataFrame(data={"Gene":colALl, "Inhibition Count": count})
    newDF.sort_values(["Gene"], inplace=True)
    # for index in newDF.index:
    #     if newDF.loc[index, 'Inhibition Count'] < 1 :
    #         newDF.drop(index, inplace=True)

    return newDF



if __name__ == '__main__':

    geneInhibitionCount = main()

    fileOut = open("GeneInhibCount.txt", "w")
    fileOut.write(geneInhibitionCount.to_string(justify="left"))
    fileOut.close()

    print("Finished Running!")
