# print(pd.options.display.max_rows)
# pd.options.display.max_rows = 9999

# for gene in totalGenesList:
#     fNameToCount = {}
#     for fIndex, fName in enumerate(dataList):
#         df = getDF(fName)
#         for dfIndex in df.index:
#             geneSubList = df.loc[dfIndex, 'Genes']
#             geneSubList = geneSubList.split(";")
#             if gene in geneSubList:
#                 if fName not in fNameToCount:
#                     count = 0
#                 else:
#                     count = fNameToCount[fName] + 1
#                 fNameToCount.update({fName: count})
#             # else: fNameToCount.update({fName: 0})
#
#     # print(fNameToCount)
#     stepDF.loc[gene] = pd.Series(fNameToCount)
# print(stepDF.info())

# print(stepDF.head().to_string(justify="left"), end="\n\n")
# stepDF.plot.barh(figsize=(10,30), stacked=True)
# plt.show()

# fig, ax = plt.subplots(figsize=(10,30))
# y_pos = np.arange(len(genesList))
# ax.barh(y_pos, numberOfPathways, align='center')
# ax.set_yticks(y_pos, labels=genesList)
# ax.invert_yaxis()  # labels read top-to-bottom
# ax.set_xlabel('Number of pathways')
# ax.set_title('Which gene affects the most pathways?')
# plt.show()

# genePathDataFrame.plot.bar(y='Number of Pathways', x='Genes')
# genePathDataFrame.plot.barh(x='Genes', figsize=(10,30), stacked=True)

# genesPathwayCorrelation = {}
# genesPathwayCorrelation.update({"Genes": totalGenesList})
# genesPathwayCorrelation.update({"Number of Pathways ": []})
# for fIndex in range(len(dataList)):
#     genesPathwayCorrelation.update({"Number of Pathways " + str(fIndex): [genesDictionary[x][1] for x in totalGenesList if genesDictionary[x][0]==fIndex]})
# genePathDataFrame = pd.DataFrame(genesPathwayCorrelation)
# genePathDataFrame.sort_values(["Number of Pathways "+str(x) for x in range(len(dataList))], inplace=True)
# genePathDataFrame.plot.barh(y='Number of Pathways', x='Genes', figsize=(10,30), stacked=True)
# plt.show()


# for gene in totalGenesList:
# fNameToCount = {}
# for fIndex, fName in enumerate(dataList):
#     df = getDF(fName)
#     # diction = {}
#     # geneslist = []
#     # numberOfPathways = []
#     # for dfIndex in df.index:
#     #     geneSubList = df.loc[dfIndex, 'Genes']
#     #     geneSubList = geneSubList.split(";")
#     #     for ig, gene in enumerate(geneSubList):
#     #         if gene not in geneslist:
#     #             geneslist.append(gene)
#
#     for dfIndex in df.index:
#         geneSubList = df.loc[dfIndex, 'Genes']
#         geneSubList = geneSubList.split(";")
#         # if gene in geneSubList:
#         for gene in geneSubList:
#             if fName not in fNameToCount: count = 0
#             else: count = fNameToCount[fName] + 1
#             fNameToCount.update({fName: count})
#             # else: fNameToCount.update({fName: 0})
#
#             # print(fNameToCount)
#             stepDF.loc[gene] = pd.Series(fNameToCount)
# print(stepDF.info())
# print(stepDF.head().to_string(justify="left"), end="\n\n")
# stepDF.plot.barh(figsize=(10,30), stacked=True)
# plt.show()


