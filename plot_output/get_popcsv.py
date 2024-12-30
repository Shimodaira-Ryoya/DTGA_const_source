import os 
import sys
sys.path.append(os.path.abspath("/Users/smdir/Documents/M1Study/DTGA"))
from output_driver.loading import load_popcsv as st

def get_pop_gene(pas,run,gen,front1=False):
    """あるrunである世代の遺伝子情報csvファイルに書かれた個体集団全ての遺伝子を取得
    Args:
        pas (str): 参照する遺伝子フォルダのパス
        run (int): 参照するrun
        gen (int): 参照する世代
    return:
        gene:二次元リスト、集団それぞれの遺伝子
    """
    df,para,gene=st.read_1gen(pas,run,gen,front1)
    return gene

def main_df_forplot(pas,runlist,genlist):
    """指定したrun,genでの遺伝子情報をすべて取得
    return para"""
    df,para=st.read_rungenlist(pas,runlist,genlist)
    return para
    
        
        

###全てのrunの最終世代の個体情報取得
