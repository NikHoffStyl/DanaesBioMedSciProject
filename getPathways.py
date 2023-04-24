import math
import os
import pandas as pd

def convertRatioToFloat(x):
    a,b = x.split("/")
    c = float(int(a)/int(b))
    return c

def getLoge(x):
    lnx = math.log(x)
    return lnx

def getDF(fileName, limPvalue=None, specifyGene=None, verbose=False):
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
        elif limPvalue is not None:
            if df.loc[index, 'P-value'] >= limPvalue:
                df.drop(index, inplace=True)
            elif df.loc[index, 'Genes'].find(specifyGene) == -1:
                df.drop(index, inplace=True)
        # elif df.loc[index, 'Term'].find(specifyGene) == -1:
        #     df.drop(index, inplace=True)

    df.drop(["Old P-value", "Old Adjusted P-value","Log P-value", "Log Adjusted P-value","P-value", "Adjusted P-value",
             "Overlap", "Odds Ratio", "Combined Score"],
            axis=1, inplace=True)

    return df



def getTableForSpecificGene(indexOfFile=None, numberOfFiles=None, combCellBioMol=False, specificGene=None, limitPval=None):
    dataList = [x for x in os.walk("data").__next__()[2] if ".txt" in x and "KEGG" not in x]
    if numberOfFiles is not None:
        if indexOfFile is not None: return -1
        else:  dataList = dataList[:numberOfFiles]
    elif indexOfFile is not None:
        if combCellBioMol:
            if "molecular" not in dataList[indexOfFile]: return -1
            else: dataList = [dataList[indexOfFile], dataList[indexOfFile].replace("molecular", "cell"),dataList[indexOfFile].replace("molecular", "bio")]
        else: dataList = [dataList[indexOfFile]]

    print(dataList)
    print("Number of Files Conbsidered: %d" % len(dataList))

    fileOut = open("geneTables.txt", "a")
    for fName in dataList:
        specificGeneDF = getDF(fName, limitPval, specificGene)
        specificGeneDF.drop(["Genes"], axis=1, inplace=True)
        if not specificGeneDF.empty:
            fileOut.write("\n\n"+"_"*100)
            fileOut.write("\nFile: %s, Gene: %s" %(fName,specificGene))
            fileOut.write("\n" + "_" * 80 + "\n")
            fileOut.write(specificGeneDF.to_string(justify="center"))

    fileOut.close()

    return True


if __name__ == '__main__':
    geneList = ["DUSP6", "MSX2", "CSRP2", "EN2", "FHL1", "IRX2", "NHLH1", "PADI3", "PAD1", "SMOC1", "CHAC14", "PENK",
                "SZL"]
    for g in geneList:
        getTableForSpecificGene(specificGene=g, limitPval=0.05)

    print("-"*10)
    print("Finished Running!")
