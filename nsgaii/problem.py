from sklearn.tree import DecisionTreeClassifier
import numpy as np
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from .tool.create_dtinfo import create_DTinfo
class DTdata:

    def __init__(self,Xtr,ytr,Xte,yte,Xname,yname,
                 mkbit,evbit,plow_mk,phigh_mk,plow_ev,phigh_ev,ev_per
                 ,secpara,cartlow,AC,SZ,MT,TS):
        #note fitness()の返り値前にあるcolumns,valuesでcsvファイルに書き込むgenedataの項目を調整する
        """dataset"""
        self.Xtr=Xtr#ndarray(np.float64)  最適化に用いる入力データセット
        self.ytr=ytr#ndarray(np.int32)    最適化に用いる出力データセット
        self.Xte=Xte#ndarray(np.float64)  精度評価に用いる入力データセット
        self.yte=yte##ndarray(np.int32)   精度評価に用いる出力データセット
        self.xn=Xname#list                入力データセットの特徴量名リスト
        self.yn=yname#list                出力データセットのクラス名リスト
        """parameter"""
        self.mkbit=mkbit    #GAで決定木作成に用いるデータサブセットを何個作るか,mkbit=3なら2^3=8個作成、g2の長さ
        self.evbit=evbit    #GAで決定木評価に用いるデータサブセットを何個作るか,evbit=3なら2^3=8個作成、g3の長さ
        self.secpara=secpara#決定木の深度情報 g4の長さ
        self.cartlow=cartlow#CARTの最大深度の最低数　
        self.plow_mk =plow_mk  #決定木作成用データサブセットに決定木作成用データセットの何割のサンプルを引っ張ってくるか
        self.phigh_mk=phigh_mk #plow~phigh
        self.plow_ev =plow_ev  #決定木評価用データサブセットに決定木評価用データセットの何割のサンプルを引っ張ってくるか
        self.phigh_ev=phigh_ev #plow~phigh
        self.ev_per=ev_per     #学習データを決定木作成用データセットと決定木評価用データセットに分けるときの後者の割合
        """fitness"""
        #遺伝子評価に用いる指標、0で不採用、1で採用
        self.AC=AC#精度：遺伝子から作成された決定木を決定木評価用データセットに通した時の正解率
        self.SZ=SZ#サイズ:作製された決定木のノード数
        self.MT=MT#遺伝子一致率
        self.TS=TS#未定
        """"""
        self.nobj=sum([self.AC,self.SZ,self.TS,self.MT])#評価関数の数
        self.xnum=self.Xtr.shape[1] #datasetの特徴量の個数
        self.genelen=self.xnum+self.mkbit+self.evbit+self.secpara#遺伝子の長さ
        
    def preprocessing(self,seed):
        """データの前処理
        Args:seed (int):シード値
        """
        self.Xmk,self.Xev,self.ymk,self.yev=train_test_split(self.Xtr, self.ytr, test_size=self.ev_per, random_state=seed)
        self.submk_index_list=create_subset_index(2**self.mkbit,self.Xmk.shape[0],self.plow_mk,self.phigh_mk,seed)
        self.subev_index_list=create_subset_index(2**self.evbit,self.Xev.shape[0],self.plow_ev,self.phigh_ev,seed)               
    
    def fitness(self, gene):
        """入力ベクトル(gene)を基に決定木モデルを作成、評価
        Args:gene (list): 入力ベクトル バイナリ列
             store_DTifo_pas:決定木構造を返すか否か、
        """
        
        """遺伝子の読み込み"""
        gene1=gene[ : self.xnum]
        if self.mkbit>=1:
            gene2=gene[self.xnum : self.xnum + self.mkbit]
        else:
            gene2=[0]
        if self.evbit>=1:
            gene3=gene[self.xnum + self.mkbit : self.xnum + self.mkbit + self.evbit]
        else:
            gene3=[0] 
        gene4=gene[self.xnum + self.mkbit + self.evbit : ]
        """遺伝子のエンコード"""
        g1_zero=judge_allzeros(gene1)#gene1がすべて0なら全て採用
        if g1_zero==True:
            gene1=[1 for _ in gene1]#gene1をすべて1に
        self.mk_i=convert_bit_dec(gene2)#採用するtr(i)のi
        self.ev_j=convert_bit_dec(gene3)#採用するte(j)のj
        depth=convert_bit_dec(gene4)     
        self.maxdepth=depth+self.cartlow #決定木の最大深度
        #使用データの洗い出し
        self.subXmk = self.Xmk[self.submk_index_list[self.mk_i],:]  
        self.subymk = self.ymk[self.submk_index_list[self.mk_i]]    
        self.subXev = self.Xev[self.subev_index_list[self.ev_j],:]  
        self.subyev = self.yev[self.subev_index_list[self.ev_j]]    
        #特徴量の選別
        self.subXmk=delete_x(self.subXmk,gene1)
        self.subXev=delete_x(self.subXev,gene1)
        #名前の選別
        self.sxn=delete_xname(self.xn,gene1)#削減された特徴量名
        #モデル生成
        self.clf=DecisionTreeClassifier(max_depth=self.maxdepth)#モデル作成，呼び出し
        self.clf.fit(self.subXmk,self.subymk)#モデル学習
        #最適化そのものには関係なし
        self.dXte=delete_x(self.Xte,gene1) #Xteデータの特徴量削減データ
        #性能評価
        self.size=self.clf.tree_.node_count#モデルのノード数
        self.ac=mtest_bycart(self.clf,self.subXev,self.subyev)#採用ev(j)での精度
        self.ac_for_te=mtest_bycart(self.clf,self.dXte,self.yte)#がちテストデータでの精度
        self.trlenp=self.subymk.shape[0]/self.ymk.shape[0]#サブセットのサンプル採用率
        self.telenp=self.subyev.shape[0]/self.yev.shape[0]
        self.trust=0
        
        if self.ac<0.8:
            self.penalty=0.8-self.ac
        else:
            self.penalty=0
    
        fitness=[]
        if self.AC==1:
            fitness.append(self.ac)
        if self.SZ==1:
            fitness.append(-self.size)
        if self.TS==1:
            fitness.append(self.trust)
        #DTinfoのための情報
        dtinfo=create_DTinfo(self.clf,self.sxn,self.yn)
                     
        
        ###遺伝子ファイルに書き込むデータ
        columns=["f1(ac)" ,"f2(size)" ,"AC3" ,"tr(i)" ,"tr_per" ,"te(j)" ,"te_per" ,"f3(trust)" ]
        values=[round(self.ac,3),self.size,round(self.ac_for_te,3),self.mk_i,round(self.trlenp,3),self.ev_j,round(self.telenp,3),round(self.trust,3)]
       
        return fitness,columns,values,dtinfo,self.MT,self.penalty
    
def split_data(X,y,te_per,seed):
    """機械学習用データを学習データとテストデータにランダムサンプリング
    Args:
        X (ndarray): 入力データ
        y (ndarray): 出力データ
        te_per (float): テストデータのサイズ
        seed (int): シード値
    Returns:
        Xtr,ytr,Xte,yte:分割されたデータ
    """
    Xtr,Xte,ytr,yte = train_test_split(X, y, test_size=te_per, random_state=seed)
    return Xtr,ytr,Xte,yte

def create_subset_index(n,d_size,plow,phigh,seed):
    """データサブセットsub(i)の要素となるデータのindexのarrayを1,2..i..nこ作る
    argument n...class:int   作成するサブセットの数
             d_size:int サブセットを作る基となるデータセットのサイズ
             plow,phigh...class:float 作成するサブセットの元の大きさに対する割合
    return   sub_index_list...class:list サブセットとして採用するサンプル番号のarrayのリスト     
    """
    np.random.seed(seed)
    sub_index_list=[]
    for i in range(n):#nのサブデータセットを作る
        s=np.random.uniform(plow,phigh)#与えられたデータセットのplow～phighの比率のサンプルを使ってサブセットを作成
        split_point=int(d_size*s)
        index=np.random.permutation(range(d_size))#配列[0,1,2,3..]を作成し要素をランダムにシャッフル
        sub_index=index[0:split_point]
        sub_index_list.append(sub_index)
    return sub_index_list

def create_BSPsubset_index(n,d_size,plow,phigh,seed):
    """bootstrappingsampling"""
    np.random.seed(seed)
    sub_index_list=[]
    for i in range(n):#nのサブデータセットを作る
        s=np.random.uniform(plow,phigh)#与えられたデータセットのplow～phighの比率のサンプルを使ってサブセットを作成
        split_point=int(d_size*s)
        index=np.random.randint(0,d_size,size=d_size)
        sub_index=index[0:split_point]
        sub_index_list.append(sub_index)
    return sub_index_list
    
def convert_bit_dec(bit_list):
    """bitのリストから変換した10進数を返す
    Args: bit_list (list,binary):0,1のリスト
    Returns: num:bit_listを10進数に直した数
    """
    num=0
    for i in range(len(bit_list)):
        num = num + 2**(len(bit_list)-i-1)*bit_list[i]
    return num
    
def delete_x(Xdata,xbit):
    """xbitにおいて0となる値があるインデックスと同じXdataのインデックスを削除
    Args:
        Xdata (ndarray): 特徴量X（入力ベクトル）のデータ、二次元配列
        xbit (_type_): 0or1のどのxiを採用するかのリスト、1で採用
    Returns:
        sXdata: 改新したデータ
    """
    adx_index=[]
    notadx_index=[]
    for i in range(len(xbit)):
        if xbit[i]==1:
            adx_index.append(i)
        elif xbit[i]==0:
            notadx_index.append(i)
    sXdata=np.delete(Xdata,notadx_index,1)#Xdataは改変されない
    return sXdata 

def delete_xname(xnames,xbit):
    """listのxnamesから,xbitが0である要素と同じインデックスを削除する
    Args:
        xnames (list):
        xbit (list):0,1のリスト
    Returns:
        sxname: 改変後xnameリスト
    """
    sxname=[]
    for i in range(len(xbit)):
        if xbit[i]==1:
            sxname.append(xnames[i])
    return sxname

def judge_allzeros(bitlist):
    """bitlistが全て0ならTrueを返す
    argument:bitlist:class:list return:TorF
    """
    allzero=not any(bitlist)
    return allzero
    
def mtest_bycart(clf,Xte,yte):
    """決定木モデルclfの性能評価  ライブラリ利用に変更2024/2/14
    Args:
        clf(decisiontreeclassfier):fit済み決定木モデル
        Xte (ndarray): テストデータの入力ベクトルX
        yte (ndarray): テストデータの出力y
    Returns:
        ac (float): 精度(テストデータの正解率)
    """
    pred = clf.predict(Xte)
    ac = accuracy_score(yte,pred)
    return ac