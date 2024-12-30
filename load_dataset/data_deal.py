import numpy as np
import category_encoders as ce

def y_encode(y,yname):
    """
    #(改良できそう)クラス名で書かれた出力データセットyをクラス名ynameの順番に対応して0,1,2...と変換する
    argument...y:class.ndarray  output dataset
               yname:class.list  class name of y
    return.....ey:class.ndarray  improvng dataset
    """
    ey=y.copy()
    for i in range(len(ey)):
        j=0
        while ey[i]!=yname[j]:
            j+=1
        ey[i]=j
    return ey

def one_hot_encode(X,xname):
    """
    #xnameで指定したカテゴリーデータである要素をonehotencodingによりエンコードする
    argument...X.class:pandas.core.frame.DataFrame ｘデータ
               xname.class:str エンコードする特徴量の名前
      return...tX.class:pandas.core.frame.DataFrame encodeしたxデータ
    """
    ce_OHE = ce.OneHotEncoder(cols=[xname])
    tX= ce_OHE.fit_transform(X)
    return tX
