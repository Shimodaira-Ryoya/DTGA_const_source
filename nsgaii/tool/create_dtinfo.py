import os
import sys
sys.path.append(os.path.abspath("/Users/smdir/Documents/M1Study/DTGA"))
from output_driver.dtinfo.dt_infomation import DTinfo

def create_DTinfo(clf,xn,yn):
    """clf情報をdtinfoに保持
    Args:
        clf (): sklearnのclfクラス
        xn (list): 使った特徴量名
        yn (list): クラス名
    Returns:
        DTi(DTinfo)
    """
    DTi=DTinfo()
    DTi.memo_from_clf(clf,xn,yn)
    return DTi

def store_DTinfo(dtinfo,pas):
    """dtinfoを保存
    Args:
        dtinfo:dtinfo
        pas:保存するフォルダとファイル名
    """
    dtinfo.fwrite_info(pas)