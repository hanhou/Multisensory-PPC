
aux= [0	50.0228	
3.2	61.8679	
6.4	74.7449	
9	76.6575	
12.8	87.8672	
25.6	99.1943	
51.2	100	
72.4	100	
76.8	100	];
coherence = aux(:,1);
coherence(1)=1;
choice2_perf = aux(:,2);
clear aux

aux =[
0	577.75	3.66956
3.2	557.763	4.26573
6.4	556.187	6.38521
9	511.15	4.64086
12.8	524.8	5.18347
25.6	405.512	1.87071
51.2	325	1.62655
72.4	304.15	0.695906
76.8	298	1.08636];
choice2_RT=aux(:,2);
clear aux

aux = [0	26.1632	
3.2	38.8985	
6.4	50.4493	
9	61.5433	
12.8	72.5176	
25.6	95.3333	
51.2	98.722	
72.4	99.6833	
76.8	99.7108];
choice4_perf = aux(:,2);
clear aux;

aux=[ 0	790.525	4.85901
3.2	739.225	7.13899
6.4	667.662	8.19986
9	677.938	6.12704
12.8	584.725	5.96848
25.6	454.95	1.68558
51.2	332.025	1.3243
72.4	316.025	0.533732
76.8	306.537	0.871783];
choice4_RT = aux(:,2);
clear aux;

subplot(221)
semilogx(coherence,choice2_perf,'o-')
hold on 
semilogx(coherence,choice4_perf,'ro-')
hold off

subplot(222)
semilogx(coherence,choice2_RT,'o-')
hold on 
semilogx(coherence,choice4_RT,'ro-')
hold off

choice2_perf=choice2_perf/100;
choice4_perf=choice4_perf/100;

save monkey_data choice2_perf choice4_perf choice2_RT choice4_RT coherence