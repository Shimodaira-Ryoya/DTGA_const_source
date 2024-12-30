import seaborn as sns
import matplotlib.pyplot as plt
import os

def plot_main_res_gen(pas,df,runlist,genlist,alpha=0.75,palette="Set1"):
    """世代ごとのプロット
    mainで生成グラフを変えたいなら変更すること
    pas...グラフ保存するフォルダまでのパス
    df...outputの解情報がまとめられたdf
    """
    
    for run in runlist:
        pltdf=df[df['gen'].isin(genlist)&df['run'].isin([run])]#dfをあるrunのあるgensに限定
        fig = plt.figure(figsize = (8,5))
        ax=sns.scatterplot(data=pltdf, x="f1(ac)", y="f2(size)", hue="gen",alpha=alpha,palette=palette)
        fig.savefig(pas+"/AC_size_scat_run"+str(run))
        plt.close(fig)
        fig2 = plt.figure(figsize = (8,5))
        ax2=sns.scatterplot(data=pltdf, x="AC3", y="f2(size)", hue="gen",alpha=alpha,palette=palette)
        fig2.savefig(pas+"/AC3_size_scat_run"+str(run))
        plt.close(fig2)
        fig3 = plt.figure(figsize = (8,5))
        ax3=sns.scatterplot(data=pltdf, x="f1(ac)", y="AC3", hue="gen",alpha=alpha,palette=palette)
        fig3.savefig(pas+"/AC_AC3_scat_run"+str(run))
        plt.close(fig3)
        
def plot_main_res_run(pas,df,runlist,genlist,alpha=0.75,palette="Set1"):
    """runごとのプロット
    mainで生成グラフを変えたいなら変更すること
    pas...グラフ保存するフォルダまでのパス
    df...outputの解情報がまとめられたdf
    """
    
    for gen in genlist:
        pltdf=df[df['run'].isin(runlist)&df['gen'].isin([gen])]#dfをあるrunのあるgensに限定
        fig = plt.figure(figsize = (8,5))
        ax=sns.scatterplot(data=pltdf, x="f1(ac)", y="f2(size)", hue="run",alpha=alpha,palette=palette)
        fig.savefig(pas+"/AC_size_scat_gen"+str(gen))
        plt.close(fig)
        fig2 = plt.figure(figsize = (8,5))
        ax2=sns.scatterplot(data=pltdf, x="AC3", y="f2(size)", hue="run",alpha=alpha,palette=palette)
        fig2.savefig(pas+"/AC3_size_scat_gen"+str(gen))
        plt.close(fig2)
        fig3 = plt.figure(figsize = (8,5))
        ax3=sns.scatterplot(data=pltdf, x="f1(ac)", y="AC3", hue="run",alpha=alpha,palette=palette)
        fig3.savefig(pas+"/AC_AC3_scat_gen"+str(gen))
        plt.close(fig3)

    
   