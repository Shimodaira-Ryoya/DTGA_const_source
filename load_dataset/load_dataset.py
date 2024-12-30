import pandas as pd
from sklearn.model_selection import train_test_split

class Fundal_Dataset:
    """学習データ管理のクラス"""
    
    def __init__(self,name,datapas="datacsv"):
        """データの名前に応じてcsvファイルを読み込み、取り込む
        Args:
            name (str): 使用するデータの名前(種類..."digits")
            datapas (str):csvファイルがあるフォルダのパス
        Return:
            self.fundal_pddata(pddataframe):読み込んだデータ
            self.X(ndarray):データのX部分
            self.y(ndarray):データのラベル部分
            self.Xname(list):特徴量の名前リスト
            self.yname(list):クラスの名前リスト
        """
        
        if name=="digits":#feature:64,class:10,instance:1797 
                          #8*8の画素で書かれた数字0-9の画像集
            file = datapas + "/digits_dataset.csv"
            self.fundal_pddata = pd.read_csv(file)
            self.Xname=self.fundal_pddata.iloc[:,1:-1].columns.tolist()
            self.yname=[0,1,2,3,4,5,6,7,8,9]
            self.X=self.fundal_pddata.iloc[:,1:-1].to_numpy()
            self.y=self.fundal_pddata.iloc[:,-1].to_numpy()
            
        if name=="drybeans":#feature16,class7,instance13611
                           #7種類の乾燥豆をカメラで撮影し，12 の寸法と 4 つの形状を得たデータ集
                           #SEKER:2025,BARBUNYA:1322,BOMBAY:522,CALI:1630,HOROZ:1928,SIRA:2636,DERMASON:3546
            file = datapas + "/drybeans_dataset.csv"
            self.fundal_pddata = pd.read_csv(file)
            self.Xname=self.fundal_pddata.iloc[:,1:-1].columns.tolist()
            self.yname=["SEKER" ,"BARBUNYA" , "BOMBAY" , "CALI" ,"HOROZ" , "SIRA" ,"DERMASON"]
            self.X=self.fundal_pddata.iloc[:,1:-1].to_numpy()
            self.y=self.fundal_pddata.iloc[:,-1].to_numpy()

    def split_data(self,trnum,seed):
        """データをトレーニングデータとテストデータにランダムに分割し、それを返す
        Args:
            trnum (int): トレーニングデータの数
            seed (int): ランダム分割のシード値
        Returns:
            Xtr,ytr,Xte,yte:トレーニングデータのXとクラス、テストデータのXとクラスのndarray
        """
        Xtr,Xte,ytr,yte = train_test_split(self.X, self.y, train_size=trnum, random_state=seed)
        return Xtr,ytr,Xte,yte 
    

        
        
        