
import pandas as pd

class parameter:
    """パラメータをcsvファイルに保存or復元を行う#要バグテスト
    """
    def __init__(self,columns):
         """パラメータをcsvファイルに保存or復元を行う
         Args:
             columns (list,str): 保存、復元するパラメータの項目名を設定
         return:
            columns,values:項目名＆値のリスト
         """
         self.columns=columns#項目名のリスト
         self.values=['None']*len(columns)#項目名に入れる値のリスト(初期値はNone)

    def fill_in_values(self,values):
        """項目名リストに対応する値をvalues(list)に入れる

        Args:
            values (list): 項目名に対応する値のリスト
        """
        if len(values)!=len(self.values):
            print("The number of values ​​doesnt match the number of item names!")
        else:
            self.values=values
                
    def store_parameter(self,pas="parameter.csv"):
        """パラメータをcsvで保存(pd.to_csvで実現)
        Args:
            pas(str): 保存ファイルパス.
        """
        mylist = list(zip(self.columns,self.values))
        df= pd.DataFrame(mylist, columns=['Columns','Values'])#データフレーム型：特徴量Columns...項目名 特徴量Values...項目に対応する値
        df.to_csv(pas)
        print("complete parameter!")
    
    def load_parameter(self,pas="parameter.csv"):
        """パラメータをcsvから読み込み
        Args:
            pas(str): 読み込むファイルパス
        """
        df = pd.read_csv(pas)
        self.columns=df.loc[:,'Columns'].tolist()
        self.values=df.loc[:,'Values'].tolist()
