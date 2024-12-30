import os
from load_dataset.load_dataset import Fundal_Dataset
from parameter_set.parameter_csv import parameter
from nsgaii.problem import DTdata
from nsgaii.nsgaii import nsgaii
from nsgaii import ea_base as ea
from plot_output import get_popcsv as getcsv, plot as pl
import timeit

def main(title='const_ac08_bit10per03090909_cartlow5', seedlist=[1,2,3,4,5], 
         dataname='digits', tr_size=1300, dataset_seed=1, 
         ngen=100, psize=100, pc=1, nvm=1, clones=False, vhigh=0.3,
         mkbit=10,evbit=10,plow_mk=0.3,phigh_mk=0.9,plow_ev=0.9,phigh_ev=0.9,ev_per=0.3,
         secpara=3,cartlow=5,
         AC=1,SZ=1,MT=0,TS=0,
         genlist=[0,20,40,60,80,100],#散布図を作る世代
         ):
    
    ###学習データの準備
    dataset=Fundal_Dataset(dataname,'datacsv')
    Xtr,ytr,Xte,yte=dataset.split_data(tr_size,dataset_seed)
    xname=dataset.Xname
    yname=dataset.yname
    
    ###最適化するもの→学習データの特徴量、サンプル、パラメータ
    dtdata=DTdata(Xtr,ytr,Xte,yte,xname,yname,
                mkbit,evbit,plow_mk,phigh_mk,plow_ev,phigh_ev,ev_per,secpara,cartlow,
                AC,SZ,MT,TS)
    ###フォルダの名前
    s='{}_({}_tr{}_seed{})_run{}'.format(title,dataname,tr_size,dataset_seed,seedlist)
    print(s)      
    ###フォルダ生成,移動
    results='output/const/' + s +'/'
    os.mkdir(results)
    os.chdir(results)#outputフォルダの中に名前ｓのフォルダを作成する、またそこにディレクトリが移動
 
    ###設定したproblem.DTmodelのパラメータをテキストに保存(再現性の確保)
    p=parameter(['dataname','tr_size','dataset_seed','mkbit','evbit','secpara','cartlow',
                'plow_mk','phigh_mk','plow_ev','phigh_ev',
                'ev_per','AC','SZ','MT','TS'])
    p.fill_in_values([dataname,tr_size,dataset_seed,mkbit,evbit,secpara,cartlow,
                      plow_mk,phigh_mk,plow_ev,phigh_ev,
                      ev_per,AC,SZ,MT,TS])
    p.store_parameter("parameter.csv")
    
    ###アルゴリズム稼働"""
    ftime = open("time.txt", "w", 1)#アルゴリズムの駆動時間を出力するファイルを作り､書き込む
    for i in seedlist:
        print("*** Run ", i, " ***")
        run = "run"+str(i)
        os.mkdir(run)
        os.chdir(run)# run i というファイルを作成しディレクトリを移動
                
        tic=timeit.default_timer()#tic-tocでアルゴリズムの駆動時間を計測
        dtdata.preprocessing(i)#データ前処理
        pop = nsgaii(evaluate = dtdata.fitness, 
           select = ea.binary_tournament_dom_cd, recombine = ea.crossover_1p, 
           mutate = ea.bit_flip_mutation, initype='binary', seed=i, 
           psize=psize, nobj=dtdata.nobj, nvar=dtdata.genelen, vlow=0, vhigh=vhigh, ngen=ngen, 
           pcx=pc, pmut=nvm/dtdata.genelen, keepclones = clones)#進化計算フェーズ
     
        toc=timeit.default_timer()#時間計測終了
        ftime.write("Nsgaii run" +str(i)+ " " + str(toc - tic) + " seconds\n")#駆動時間の書き込み
        os.chdir("..")#ディレクトリを一つ前に戻す
        
    ftime.close()#駆動時間に関するファイルを閉じる
    
    
    ###グラフ生成
    os.mkdir("graph")
    os.mkdir("graph/gen_graph")
    os.mkdir("graph/run_graph")
    para=getcsv.main_df_forplot(".",seedlist,genlist)
    pl.plot_main_res_gen("graph/gen_graph",para,seedlist,genlist)
    pl.plot_main_res_run("graph/run_graph",para,seedlist,[100])

if __name__ == '__main__':
    main()