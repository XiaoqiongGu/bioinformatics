@author: Xiaoqiong Gu
@update time: 2022 Apr

### commonly used headings

	%matplotlib inline
	import numpy as np
	import pandas as pd
	import seaborn as sns;sns.set(color_codes=True)
	import matplotlib.pyplot as plt
	from scipy import stats
	from scipy.stats import zscore
	
	import glob
	import os
	import argparse
	
	# dendrogram importing
	from matplotlib import pyplot as plt
	from scipy.cluster import hierarchy

## pandas commonly used function/scripts
### File import, export

import files from excel or csv or txt

	x = pd.read_excel("x.xlsx",sheet_name='x')
	x = pd.read_csv('x.txt',sep='\t',header=None/0,index_col=False/0/'columnA') 
	#no row as header, #no column/the 1st column/column A as index header

export dataframe to csv or txt

	df.to_csv("df.csv",index=False)
	df.to_csv("df.txt",sep="\t", index=False) 

view all the column content

	pd.set_option('display.max_columns', None)

### file processing

[merge two files based on share column](https://stackoverflow.com/questions/40468069/merge-two-dataframes-by-index/40468090#40468090)

	pd.merge(df1, df2, left_index=True, right_index=True, how='inner') 
	pd.merge(df1, df2, left_on='df1_key', right_on='df2_key',how="outer")
	how type: inner, outer, right, left
	merge by two columns: left_on=['A1','A2'], right_on = ['B1','B2'])

melt the df, wide to long format needed for seaborn plot in hue setting

	pd.melt(X.reset_index(),id_vars='p',value_vars = ['S1','S2' ..'S30'],var_name="subject",value_name='rel')  
	frame : 你想要更動的 DataFrame。
	id_vars: 可使用 tuple、list、或 ndarray，用以設定不想要被轉換的欄位。
	value_vars: 可使用 tuple、list、或 ndarray，用以設定想要被拆解的欄位。 如果省略則拆解全部欄位。
	var_name : 轉換後 id 的名稱。如果省略則設定為原本 DataFrame 的欄位名稱或是 variable。
	value_name : 轉換後 value 欄位的名稱。如果省略則顯示原本 DataFrame 的欄位名稱或 value。
	col_level : 可使用 int、string。如果 columns 是 MultiIndex，則使用該參數來進行選擇。


pivot the melted dataframe, long to wide format
	
	df.pivot(index=None, columns=None, values=None)
	df.pivot(index='arg', columns='days', values='rel')


create a new column based on one column in the dataframe, apply lambda function

	df['new column name'] = df['column name'].apply(lambda x: 'value if condition is met' if x condition else 'value if condition is not met')
	df['days'] = df['subject'].apply(lambda x: int(x.split('D')[0]))
	dfmm['aro_accession'] = dfmm['arg'].apply(lambda x: x.split('|')[2])
	df['MajorAllele'] = df.apply(lambda x: x.ALT if x.SNP_frac >=0.5 else x.REF,axis=1) # note havec to put axis = 1

apply certain functioin on series
	
	df['email'].apply(len)

	def update_email(email):
		return email.upper()

	df['email'].apply(update_email) 
	# just want passing the function instead of excutive the function, so do not need to add the ()

	lamda function -> if the function is not too complicated  
	df['email'].apply(lambda x: x.lower())

create a new column based on two columns in the dataframe, apply lambda function

	def f(x):    
   		return x[0] + x[1] 
    df.apply(f, axis=1) #passes a Series object, row-wise


	def f(sta,end):
    	return mylist[sta:end+1]
	df['col_3'] = df.apply(lambda x: f(x.col_1, x.col_2), axis=1)

apply on the dataframe
	
	df.apply(len, axis='columns')
	https://datatofish.com/if-condition-in-pandas-dataframe/


replace nan with 0

	df.fillna(0)
	df.replace(np.nan,0)


get the relative abundance of the absolute OTU table

	bacteria_rel = bacteria.divide(bacteria.sum(axis = 0),axis =1)
	otu97.groupby(['taxonomy']).sum()
	axis = 0 --> based on rows 
	axis = 1 --> based on columns

	DataFrame.groupby function on columns, set axis=1 (iterate over each column), apply mean function for each group.
	df.groupby(by=df.columns, axis=1).mean()


get the mean or average of the dataframe

	df.mean(axis=0)



### file format conversion

convert set to list

	list(set)

convert series to df

	y = pd.Series.to_frame(x)

convert index to list
	
	pd.index.to_list()


convert the column string to interger in pandas dataframe
	
	df['df Column'] = df['df Column'].astype(int)
	df['df Column'] = pd.to_numeric(df['df Column'])

Convert columns to string in Pandas

	total_rows['ColumnID'] = total_rows['ColumnID'].astype(str)

#### selecting columns

select columns based on columns names containing a specific string in pandas

	df.loc[:, df.columns.str.startswith('alp')]
	df.filter(like='keyword').head
	df.filter(regex='keyword').head

select multiple columns based on the columns string that contain keyword D0

	df.filter(regex='D0')
	df.filter(regex='S1D|S2D')

drop certain columns in DataFrame

	df.drop(df.columns[[1, 2]], axis=1, inplace=True)

drop by name

	df.drop(['B', 'C'], axis=1, inplace=True) # axis=1, column; axis=0, row
	df.drop(columns=['B'], inplace=True) 

drop duplicates
	
	df.drop_duplicates()

select duplicates
	
	df.duplicates()
	
select single column

	df['X'].head() -> return the tuple object
	df[['X']] -> return the dataframe, 2-dimensional data
	df.X.head()

[select single/multiple columns](https://www.coder.work/article/85145)

1. loc 方法
	`df_new = df.loc[:, 'a']`
	`df_new = df.iloc[:,0:2]`  #we need the : to select the colomn index

2. 单层中括号
	`df_new = df['a']`
	`df[['a','b']]`

3. 最简单的方法
	`df_new = df.col1`

#### selecting rows

select rows based on row names in pandas

	df.loc[["y"]]

> Note, two square bracket instead of one square bracket

split the columns

	df['V'].str.split('-',expand=True)
	df['v'].str.split('-').str[-1] => only keep the last column 
	df[['V1','V2']]= df['V'].str.split -> return the series, once expand, can return the dataframe, if not expand, perhaps return the list? 


sort values and choose the top N index
	pd.sort_values(by='A',ascending=False)[:N].index()

select columns containing partial string
	df[df['A'].str.contains("keyword")]

### index

1. you can set column A and B simutaneously as the index column
	df.set_index(['A','B'],inplace=True)
2. return the index columns A and B as the normal columns A and B
	df.reset_index(inplace=True) 
3. sort the index
	df.sort_index(ascending=True/False)
	df.loc[index,'columnA'] # why using the index is useful? 

### filter
1. df how to filter based on column value (boolean variables with True and False in it)

		df[df["a"]==0]
		df.loc[df["a"]==0] # you can use .loc to grab the specific column or make a new column

		filt = (df["a"]==0)
		df.loc[filt, 'columnA']

2. &,|, ~  -> operator 
		
		filt = (df["a"]==0) & (df['b'] ==0)
		df.loc[\~filt, 'columnA] # use tilt to negate the filter
3. df['country'].isin(countries)
4. str.contains('keyword',na=False) #we will not do anything with those nan


### update rows and columns
#### update column names

1. rename columns

		df.rename({'A':'A1'},axis = columns, inplace=True) #you can specify certain column names
2. use the list comprehension: to add a prefix or suffix to each column name, 
		
		df.columns = [str(col) + '_x' for col in df.columns]
		df.columns = [x.upper() for x in df.columns] 
		df = df.add_suffix('_some_suffix')
3. str.replace
	
		df.columns.str.replace(' ','_')
3. to remove a prefix or suffix to each column name
	
		import pandas as pd
		def drop_prefix(self, prefix):
			self.columns = self.columns.str.lstrip(prefix)
			return self
		Then you can use it as with inverse method already implemented in pandas add_prefix:
		pd.drop_prefix('myprefix_')

#### update rows in the content
1. df.loc[2, ['A','B']] 
2. update all the rows: 
	
	df['email'] = df['email'].str.lower()


Grouping and Aggregating

	alue_counts(normalize=True) # very useful and maybe can be used to normalize the OTU features

split, apply, combination  


rearrange the row order or column order

	#the row order
	df.reindex()
	#the column order
	cols = df.columns.tolist() -> get the list of columns and then rearrange
	df = df[cols]

reinterate the list and select the specific locations in the list

	list1= []
	list2 = [x[0] for x in list1]


join, which is left join by default:

	df1.join(df2)

combine two series 	or dataframe

concat, which is outer join by default:

	pd.concat([df1, df2], axis='rows') 

combine two dataframe by joint rows

	pd.concat((df1,df2),axis='columns',sort=False)


### numpy
Returns the indices that would sort an array.
	
	A.argsort()  # increasing order
	(-A).argsort() # decreasing order

rearange array according to the index order
	
	A[idx]

### useful function

Pandas Count Occurrences in Column – i.e. Unique Values

	df['column'].value_counts()
	df['column'].value_counts(dropna=False) #including missing values
	df['column'].value_counts().ColValue #specific value
	df.groupby('column').count()

Get value counts for multiple columns at once in Pandas DataFrame

	subdf[['qstart','qend']].apply(pd.Series.value_counts)


### [regular expression usage]('https://www.regular-expressions.info/index.html')
	xx = ''
	r1 = re.findall(r'VL\w*$',xx)



## adonis VS PERMANOVA
http://qiime.org/scripts/compare_categories.html



## Seaborn and matplotlib

	import seaborn as sns; sns.set()
	import matplotlib.pyplot as plt
	sns.set_style(style="white")
	sns.set(font_scale=2)
	
set the figure size

	fig, ax = plt.subplots(figsize=fig_dims)


make subplots 

	fig, ax = plt.subplot(nrows=n, ncols=m, plotNum)
	

iterate over many subplots
	
	fig, axes = plt.subplots(nrows=n, ncols=m)
	axes is then an array of n rows and m columns. One could then iterate over the rows and columns independently, or over the flattened version of the axes array.

	x = np.arange(11)
	y = np.random.rand(len(x), 9)*10

	fig, axes = plt.subplots(3,3, sharex=True, sharey=True)
	for i, ax in enumerate(axes.flatten()):
    	ax.bar(x, y[:,i], color=plt.cm.Paired(i/10.))
	plt.show()

	from matplotlib.ticker import PercentFormatter
	for i,ax in enumerate(axes.flatten()):
    	ax.bar(x,y.iloc[i,:],color=colorcodes[i],width=0.5)
    	ax.set_ylabel(y.iloc[i,:].name.split(';')[-1])
    	ax.yaxis.set_major_formatter(PercentFormatter(1,decimals=0)) 

gridspec

	https://matplotlib.org/stable/tutorials/intermediate/gridspec.html


## seaborn make the picture - do a function 一劳永逸
	def create_heatmap(dataframe, title):
	import seaborn as sns
	import matplotlib.pyplot as pyplot	
	sns_plot=sns.clustermap(dataframe,mask=False,row_cluster=True,linewidth=.05,linecolor="grey",col_cluster=False,standard_scale=0,cmap="Greys",square=True,figsize=(12,(.2*(len(dataframe.index)))))
	sns_plot.savefig("{}.png".format(title))
	return sns_plot


## color palette

	colorbrew2
	 
extract RGB or 6 digit code from seaborn palette

	print(sns.color_palette('Set3').as_hex())
	sns.color_palette('Set2')

Set your custom color palette
	
	colorcodes = ['#8dd3c7', '#fdb462', '#bebada', '#fb8072', '#80b1d3','#ffed6f',  '#b3de69', '#fccde5', '#d9d9d9', '#bc80bd']
	sns.set_palette(sns.color_palette(colorcodes))
	sns.color_palette()

set a darker version of the custom color palette

	def lighten_color(color, amount):
    	"""
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

	darken_colorcodes = [lighten_color(x,2) for x in colorcodes]
	sns.color_palette(darken_colorcodes)
	https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib/49601444

	darken_colorcodes = ['#368c7d','#be6502','#7d75b5','#d61b06','#2b5b7d','#dec200','#54761a', '#f99bcb', '#b3b3b3', '#512a52'])

my fav color palette
* *red + blue* 
	
		palette = ['#e41a1c','#377eb8'] 
		palette = ['#fb8072', '#80b1d3'] # light version
* *shades of gray* [sources](https://www.w3schools.com/colors/colors_shades.asp)

		'#383838'

* *stacked bar plots with many species*

		colorcodes = ['#EF926C','#7CC4A7','#90A4CA','#DA90C4','#B2DB68','#6363A7','#E2C799','#84CBEB','#DF6C67','#878428','#F9E053','#6C3C17','#C6C6C6','#848484']

## log scale in scatter plot
	plt.xscale('log')
	plt.yscale('log')
	

### matplotlib.pyplot
	import matplotlib.pyplot as plt
	plt.style.use('') # print(plt.style.available)
	%matplotlib inline
>	在使用jupyter notebook 的时候，经常会用到%matplotlib inline

	其作用就是在你调用plot()进行画图或者直接输入Figure的实例对象的时候，会自动的显示并把figure嵌入到console中。
	注意的地方
	但有一个不太方便的地方，当你调用fig1 = plt.figure(1); 如果再次调用plt.figure(1)的时候会产生新的Figure实例对象，而且每次plt.gca()或者plt.gcf()都会产生不同的对象。
	如果不加上%maplotlib inline 的话，每次figure的显示都需要plt.show();

histogram plot

	(n, bins, patches) = plt.hist(x, bins=45, label='hst')



swarmplot
	
	sns.swarmplot(x='Group',y='value', data=X,color=".25",ax=ax0, size=10)

boxplot
	
	sns.boxplot(x='Group',y='value', whis=1.5, data=X,ax=ax0, palette = ['#e41a1c','#e41a1c','#377eb8','#377eb8']) # showfliers=False determined whether you wanna show outliers or not

lineplot

	sns.lineplot(data=otu97_rel_p_AAD_avg.T, ax=ax0,linewidth=3)

areaplot

	df.plot(kind='area',figsize=(5,5),linewidth=0,colormap='Set3',ax=ax_loc)

barplot

option 1	

	https://matplotlib.org/3.1.1/gallery/lines_bars_and_markers/barchart.html#sphx-glr-gallery-lines-bars-and-markers-barchart-py
	grouped bar plot

	x = np.arange(len(labels))  # the label locations
	width = 0.35  # the width of the bars

	fig, ax = plt.subplots()
	rects1 = ax.bar(x - width/2, men_means, width, label='Men', color=colorcodes[i],yerr=error)
	rects2 = ax.bar(x + width/2, women_means, width, label='Women')

option 2

	df.plot.bar(stacked=True, linewidth=1, width=0.8,edgecolor='black',ax=ax)

	# insert labels above every columns
	rects = ax.patches
	# Make some labels.
	labels = ['4033','227','240','266','278','2726','421','412','283','314','387']
	for rect, label in zip(rects, labels):
    	#height = rect.get_height()
    	ax.text(rect.get_x()+rect.get_width()/2, 1.05,label, ha="center", va="bottom") # this is assuming it is the stacked bar plots when the heights are fixed, otherwise can use height variable



scatterplot

	sns.scatterplot(data=metadata,x='total_number_of_episodes',y='maximum_BS',hue='Group',s=size,palette=[ '#377eb8','#e41a1c'])

clusterplot
	
	sns.clustermap(corr_df,row_cluster=False, col_cluster=False,mask=mask ,cmap="RdBu_r")
	corr_df=.T.corr(method='spearman')

	ax = sns.clustermap(corr_df,cmap="RdBu_r") #row_cluster=False, col_cluster=False,mask=mask,
	corr_df.index[ax.dendrogram_row.reordered_ind]

pie plot
	
	plt.pie(x=otu97_rel_g_rumino_avg.loc[top5_rumino_g]['Ruminococcaceae'], labels=labels, autopct="%.1f%%", explode=[0.03]xN,textprops={'fontsize': 20}, pctdistance=0.5) # N equals to the length of x


heatmap plotting

	import numpy as np; np.random.seed(0)
	import seaborn as sns; sns.set()
	uniform_data = np.random.rand(10, 12)
	ax = sns.heatmap(uniform_data)

### rc parameters

	font = {'family': 'arial',
        'weight': 'light',
        'size': 16}
	plt.rc('font', **font)
	plt.rcParams["axes.edgecolor"] = "black"
	plt.rcParams["axes.linewidth"] = 1.8   ## setting the frame width

### legend
	plt.legend(fontsize=26,markerscale=2)
	plt.plot(faecaliqpcr_aad_normalized,color='#fb8072', linewidth=5,linestyle=':',label='$\it{Faecalibacterium}$ $\it{Prausnitzii}$-qPCR')
	$\it{Faecalibacterium}$ $\it{Prausnitzii}$
	
	$\it{}$

	ax.legend().set_title('')
	ax.get_legend().remove()
	ax.set_title('AAD group (N=13)')
	ax.legend(loc=5, bbox_to_anchor=(1.7, 0.5),fontsize=20,markerscale=2)

handels and labels for seaborn: https://stackoverflow.com/questions/51579215/remove-seaborn-lineplot-legend-title

	handles, labels = ax.get_legend_handles_labels()
	ax.legend(handles=handles[1:], labels=labels[1:])

[handels location/legend location](https://stackoverflow.com/questions/44413020/how-to-specify-legend-position-in-matplotlib-in-graph-coordinates)

	handles, labels = ax[0].get_legend_handles_labels()
	fig.legend(handles, labels, bbox_to_anchor=(0.4,0),loc='center')

	h,l = ax.get_legend_handles_labels()
	plt.legend(h[0:3],l[0:3],bbox_to_anchor=(1.05, 1), loc='upper left, borderaxespad=0., fontsize=26,markerscale=2)
	https://stackoverflow.com/questions/4700614/how-to-put-the-legend-out-of-the-plot

plot x horizontal line

	plt.axhline(y=N, color='r', linestyle='--')

plot y horizontal line

	plt.ayhline(y=N, color='r', linestyle='--')

plt tick rotation

	plt.xticks(rotation = 90)

### ax parameter setting
change x ticklabel

	for tick_labels in ax.axes.get_xticklabels():
    	tick_text = tick_labels.get_text()
	
	for tick_labels in ax.xaxis.get_ticklabels():
    	print(tick_labels)


	plt.draw()
	tick_text = [t.get_text().split('|')[-1]  for t in ax.get_xticklabels()]
	ax.set_xticklabels(tick_text)

adding the grid line

	ax.grid(axis='both', which='major', ls='--', alpha=0.4)

set the ax spine linewidth

	for axis in ['top','bottom','left','right']:
  		ax.spines[axis].set_linewidth(4)


ax parameters

	ax1.set_xlabel('Baseline samples',fontsize=16)
	ax1.set_xlabel('') #set an empty string and to remove the xlabel
	ax1.set_ylabel(r'$\alpha$ Diversity'+'\n(chao1)',fontsize=16)
	ax1.tick_params(axis='x',labelsize=14,direction="in")
	ax1.tick_params(axis='y',labelsize=14,direction="in")
	ax.tick_params(axis='both',labelsize=16,direction="in",length=x,width=y)
	ax.yaxis.tick_left()
	ax.xaxis.tick_bottom()
	
	ax.set_xticks(np.arange(len(xtickid))) # this is to set the tick so that all the labels can be shown in the set_xticklabels
	ax.set_xticklabels(['0','1','2','3','7','14','28']) #reset xtick labels
	ax.set_xlim([xmin,xmax])
	ax.set_ylim([ymin,ymax])

Y-axis ticks on right side of plot

	ax.yaxis.tick_right()
	add_stat_annotation(ax, data=subdf, x='Day', y='Abundance',line_height=0.01,hue='Group',order=day_order,box_pairs = box_pairs, test='Mann-Whitney', comparisons_correction='bonferroni' or 'None', text_format='star', loc='outside', verbose=2, fontsize=20)


[Show decimal places and scientific notation on the axis of a matplotlib plot](https://stackoverflow.com/questions/25750170/show-decimal-places-and-scientific-notation-on-the-axis-of-a-matplotlib-plot)
	import matplotlib.ticker as mticker
	f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
	g = lambda x,pos : "${}$".format(f._formatSciNotation('%1.10e' % x))
	ax1.get_yaxis().set_major_formatter(mticker.FuncFormatter(g))

### Hide non-significant annotations
https://github.com/webermarcolivier/statannot/issues/50


### cmap database

	supported values are 'Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'crest', 'crest_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'flare', 'flare_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'hot', 'hot_r', 'hsv', 'hsv_r', 'icefire', 'icefire_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'mako', 'mako_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'rocket', 'rocket_r', 'seismic', 'seismic_r', 'spring', 'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'turbo', 'turbo_r', 'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'viridis', 'viridis_r', 'vlag', 'vlag_r', 'winter', 'winter_r'


# statistics
normalization

	from scipy.stats import zscore
	df.apply(zscore)

describe the dataframe

	df.describe()


choose different types

 	t-test_ind, t-test_welch, t-test_paired, Mann-Whitney, Mann-Whitney-gt, Mann-Whitney-ls, Levene, Wilcoxon, Kruskal


	mwu statistics is for the two independent samples
	Wilcoxon signed-rank test is for the paried/matched samples
	https://stats.stackexchange.com/questions/113936/what-is-the-difference-between-the-mann-whitney-and-wilcoxon-rank-sumtest

[kstest](https://www.statology.org/kolmogorov-smirnov-test-python/)
	
	The Kolmogorov-Smirnov test is used to test whether or not or not a sample comes from a certain distribution.
	from scipy.stats import kstest   # one sample
	from scipy.stats import ks_2samp # two samples
	pvalue = []
		for i in range(df.shape[0]):
    		pvalue.append(ks_2samp(df1.iloc[i, :], df2.iloc[i,:]).pvalue)

[Binomial distribution](https://www.statology.org/binomial-test-python/)
	
	stats.binom_test(samp_x, n=samp_n,p=hypo_p, alternative='less')


[multiple hypothesis test](https://towardsdatascience.com/multiple-hypothesis-testing-correction-for-data-scientist-46d3a3d1611d)

	from statsmodels.stats import multitest 
	mw_stats_fdrbh=multitest.multipletests(mw_stats[0], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
	mw_stats_fdrbh[0]
	sig_fdr=np.where(mw_stats_fdrbh[0])
	mw_stats_fdr_bh[1][sig_fdr]

### python regex
https://www.programiz.com/python-programming/regex

### python 3
	string formatting
	name ='giacomo'
	number = 4.3
	print('%s %s %d %f %g' % (name, number, number, number, number))





### datetime module
convert datetime to date

	datetime.datetime.now().date()


[the inverse of datetime.isocalendar()](https://stackoverflow.com/questions/304256/whats-the-best-way-to-find-the-inverse-of-datetime-isocalendar)

	Python 3.8 added the fromisocalendar() method:
	>>> datetime.fromisocalendar(2011, 22, 1)
	datetime.datetime(2011, 5, 30, 0, 0)
	
	>>> datetime.strptime('2011 22 1', '%G %V %u')
	datetime.datetime(2011, 5, 30, 0, 0)



### Others

secondary axis

	ax = sns.lineplot(data=pdd, x="pos_adj", y="total_cov",color='b',lw=2)
	ax2 = ax.twinx()
	sns.barplot(data=pdd,x='pos_adj',y='fraction_break',color='g',ax=ax2)

run bash scripts in python scripts and get the output of bash scripts

	process = subprocess.Popen([ 'cat '+ sys.argv[1]+'|cut -f2|sort -u|wc -l'], shell=True, stdout=subprocess.PIPE)
	ecoli_c0 = int(process.communicate()[0])
	process.communicate() return stdout and stderr


### if __name__=='main'
	def main():
		{scripts}
	if __name__=='main':
		main()

	__name__ 是当前模块名，当模块被直接运行时模块名为 __main__ 。这句话的意思就是，当模块被直接运行时，以下代码块将被运行，当模块是被导入时，代码块不被运行。



### color hex code
color choice

	blue: #A6BAE4, #A7BAE2
	light blue: #DADFEB
	dark blue: #D0D4E5
	
	green: #ADBD93, #A5B48B
	Left green: #E8F0E7
	Right Green: #E5EFF1,#8DA1A6
	
	orange: #f2c98a, #D9B47B, #DEB87F
	red: #9A7484
	Light gray: #E7E9E9
	Dark green: #96AD9F, #789485
	Gray: #6D6E71
