# 教程来源 
# https://projects.volkamerlab.org/teachopencadd/talktorials/T001_query_chembl.html#

# 加载必要的package
import math
from pathlib import Path
from zipfile import ZipFile
from tempfile import TemporaryDirectory
import numpy as np
import pandas as pd
from rdkit.Chem import PandasTools
from chembl_webresource_client.new_client import new_client
from chembl_webresource_client.settings import Settings
Settings.Instance().TIMEOUT = 10
from tqdm.auto import tqdm

import pubchempy as pcp


# 初始化活性的api
bioactivities_api = new_client.activity

# 本方法需要化合物的CHEMBL ID因此，可以通过inchiKey进行匹配
# inchiKey可以通过PubChem Identifier Exchange Service基于CID获取


jiandan = False
if jiandan:
    #简单的方法，直接有现成的inchiKey
    inchi_key = pd.read_csv("CHEMBL_ID.txt",sep="\t")
else:
    print("没有那么简单的事情~")
    # 【提前准备好】读取chembl全部注释的化合物,可以从官网获取
    CHEMBL_chemical = pd.read_csv("chembl_30_chemreps.txt.gz",sep="\t")
    # 【提前准备好】读取含有名字和pubchemCID的初始文件,注意列名需要一致
    pubchem_cid = pd.read_csv("pubchem_cid.txt",sep="\t")
    
    
    df1 = pcp.get_properties(['inchikey'], list(pubchem_cid['Pubchem Cid']) , as_dataframe=True)
    df2 = df1.stack().reset_index().drop("level_1",axis=1)
    df2.columns = ['Pubchem Cid', 'standard_inchi_key']
    df3 = pd.merge(pubchem_cid,df2,how="left")
    inchi_key = pd.merge(df3,CHEMBL_chemical,how="left")
    # inchi_key.to_csv("CHEMBL_ID.txt",sep="\t",index=False)

# 构建请求api，准备下一步请求
bioactivities = bioactivities_api.filter(
    # target_chembl_id=chembl_id,  # 搜索靶点
    molecule_chembl_id__in=list(inchi_key["chembl_id"].dropna()),
    # molecule_pref_name = "Baicalin", # 基于名字搜索，但是需要全部大写
    type="IC50", 
    # relation="=", #指标参数
    # assay_type="F" # chembl官网有说明，靶点一般是B，IC50一般是F
).only( # only是简化功能，即我只需要下面的信息，可以加快运行速度
    "activity_id",
    "assay_chembl_id",
    "assay_type",
    "assay_description",
    "assay_type",
    "bao_label",
    "molecule_chembl_id",
    "molecule_pref_name",
    "type",
    "standard_units",
    "relation",
    "standard_value",
    "target_chembl_id",
    "target_organism",
    "target_pref_name"
)
# 看一下数据量有多少
print(f"Length and type of bioactivities object: {len(bioactivities)}, {type(bioactivities)}")


# 表格化数据,这一步会联网，视数据量的多少时间不确定，多的话时间会比较久
# 这个democode就很久
bioactivities_df = pd.DataFrame.from_records(bioactivities)

# 查看换一下列名，方便下一步提取
bioactivities_df.columns

bioactivities_df_all = pd.merge(inchi_key,bioactivities_df,left_on='chembl_id',right_on='molecule_chembl_id')
bioactivities_df_all.to_excel("all.xlsx")
# 数据提取，提取想要的细胞系啥的,这一步可以自定义
bioactivities_df1 = bioactivities_df[(bioactivities_df['bao_label'] == 'cell-based format') 
                                     & (bioactivities_df['relation'] == '=') 
                                     & (bioactivities_df['target_organism'] =='Homo sapiens')
                                     & (bioactivities_df['target_pref_name'] =='PC-3')
                                    ]
                                    
# 结果合并并输出
res = pd.merge(inchi_key,bioactivities_df1,left_on='chembl_id',right_on='molecule_chembl_id')


# 输出结果                                    
res.to_excel("res.xlsx")



# print(f"DataFrame shape: {bioactivities_df.shape}")
# bioactivities_df.head()
# 
# print(f"Length and type of bioactivities object: {len(bioactivities)}, {type(bioactivities)}")
# bioactivities_df["units"].unique()


###### 获取 target 或者是 molecule的结构
# targets_api = new_client.target
# compounds_api = new_client.molecule
# uniprot_id = "P00533"
# # Get target information from ChEMBL but restrict it to specified values only
# targets = targets_api.get(target_components__accession=uniprot_id).only(
#     "target_chembl_id", "organism", "pref_name", "target_type"
# )
# print(f'The type of the targets is "{type(targets)}"')
# 
# 
# targets = pd.DataFrame.from_records(targets)
# targets
# target = targets.iloc[0]
# target
# 
# chembl_id = target.target_chembl_id
# print(f"The target ChEMBL ID is {chembl_id}")
# # NBVAL_CHECK_OUTPUT
# 
# compounds_provider = compounds_api.filter(
#     molecule_chembl_id__in=list(inchi_key["chembl_id"].dropna())
# ).only("molecule_chembl_id", "molecule_structures")
# 
# compounds = list(tqdm(compounds_provider))


