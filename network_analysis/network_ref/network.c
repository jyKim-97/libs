#include <limits.h>


void SF_configuration(int* linkn, int n_node, double avg_link, double gamma)	{

int i,j;

for(i=0;i<n_node;i++)	{
	linkn[i]=0;
}

for(i=0;i<n_node;i++)	{
	linkn[i]=(int)((double)pow((16./(1.-ran4())),0.5)-3.);
}

}

void SF_static(int** edge, int* linkn, int n_node, double avg_link, double gamma)
{

int i,j;
int *prob_pnt,*poor,*rich;
int npoor,nrich,tot_degree,current_degree,des,dep,tmp,int_tmp;
double *prob;
double alpha,sum_prob,prob_critical,db_tmp;

alpha=1.0/(gamma-1.0);

for(i=0;i<n_node;i++)
{
	edge[i]=(int*)realloc(edge[i],sizeof(int)*1);
	edge[i][0]=(-1);
	linkn[i]=0;
}

prob_pnt=(int*)malloc(sizeof(int)*n_node);
prob=(double*)malloc(sizeof(double)*n_node);


rich=(int*)malloc(sizeof(int)*1);
poor=(int*)malloc(sizeof(int)*1);
rich[0]=poor[0]=0;

for(i=0,sum_prob=0.;i<n_node;i++)
{
	prob[i]=pow((double)(i+1),-alpha);
	prob_pnt[i]=i;
	sum_prob+=prob[i];
}

for(i=0;i<n_node;i++)	{prob[i]=prob[i]/(double)sum_prob*(double)n_node;}

prob_critical=1.0;

for(i=0,npoor=nrich=0;i<n_node;i++)
{
	if(prob[i]>prob_critical)
	{
		nrich++;
		rich=(int*)realloc(rich,sizeof(int)*nrich);
		rich[nrich-1]=i;
	}
	else if(prob[i]<=prob_critical)
	{
		npoor++;
		poor=(int*)realloc(poor,sizeof(int)*npoor);
		poor[npoor-1]=i;
	}
}

while(1)
{
	prob[rich[nrich-1]]-=prob_critical-prob[poor[npoor-1]];
	prob_pnt[poor[npoor-1]]=rich[nrich-1];
	npoor--;
	poor=(int*)realloc(poor,sizeof(int)*npoor);
	if(prob[rich[nrich-1]]<prob_critical)
	{
		npoor++;
		poor=(int*)realloc(poor,sizeof(int)*npoor);
		poor[npoor-1]=rich[nrich-1];
		nrich--;
		rich=(int*)realloc(rich,sizeof(int)*nrich);
	}

	else if(prob[rich[nrich-1]]==prob_critical)
	{
		nrich--;
		rich=(int*)realloc(rich,sizeof(int)*nrich);
	}

	if(npoor<1 || nrich<1)	{break;}
}

tot_degree=(int)(((double)avg_link/2.0)*n_node);
current_degree=0;

while(1)
{
	tmp=(int)(ran4()*n_node);
	if(prob[tmp]>ran4())	{dep=tmp;}
	else {dep=prob_pnt[tmp];}


	tmp=(int)(ran4()*n_node);
	if(prob[tmp]>ran4())	{des=tmp;}
	else {des=prob_pnt[tmp];}

	if(dep!=des)
	{
		for(tmp=0,i=0;i<linkn[dep];i++)	{if(des==edge[dep][i])	{tmp=1;break;}}

		if(tmp==0)
		{
			linkn[dep]++;linkn[des]++;

			edge[dep]=(int*)realloc(edge[dep],sizeof(int)*linkn[dep]);
			edge[des]=(int*)realloc(edge[des],sizeof(int)*linkn[des]);

			edge[dep][linkn[dep]-1]=des;
			edge[des][linkn[des]-1]=dep;

			current_degree++;
		}
	}
	if(current_degree==tot_degree)	{break;}
}

free(poor);free(rich);
free(prob);free(prob_pnt);


}



void Random_ErdosRenyi(int** edge, int* linkn, int n_node, double avg_link)
{

int i,j;
int tot_degree,current_degree,dep,des,tmp;

for(i=0;i<n_node;i++)
{
	edge[i]=(int*)realloc(edge[i],sizeof(int)*1);
	edge[i][0]=(-1);
	linkn[i]=0;
}


tot_degree=(int)(((double)avg_link/2.0)*n_node);

current_degree=0;

while(1)
{
   dep=(int)(ran4()*n_node);
   des=(int)(ran4()*n_node);

   if(dep!=des)
   {
      for(tmp=0,i=0;i<linkn[dep];i++)  {if(des==edge[dep][i])  {tmp=1;break;}}

      if(tmp==0)
      {
         linkn[dep]++;linkn[des]++;

         edge[dep]=(int*)realloc(edge[dep],sizeof(int)*linkn[dep]);
         edge[des]=(int*)realloc(edge[des],sizeof(int)*linkn[des]);

         edge[dep][linkn[dep]-1]=des;
         edge[des][linkn[des]-1]=dep;

         current_degree++;
      }
   }
   if(current_degree==tot_degree)   {break;}
}

}

int Burning_algorithm(int** edge, int* linkn, int* flags, int n_node)
{

int i,j,tmp;
int *nn,*sizecl;
int r_node,cluster,largest_cluster;

cluster=(-1);
nn=(int*)malloc(sizeof(int)*n_node);
sizecl=(int*)malloc(sizeof(int)*(n_node+1));
r_node=0;


for(i=0;i<n_node;i++)   {flags[i]=nn[i]=(-1);}
for(i=0;i<n_node+1;i++)	{sizecl[i]=0;}

for(i=0;i<n_node;i++)
{
   if(flags[i]==(-1))
   {
      cluster++;
      flags[i]=cluster;
      sizecl[cluster]=1;
      tmp=0;
      nn[tmp]=i;

      while(1)
      {
         for(j=0;j<linkn[nn[0]];j++)
         {
            if(flags[edge[nn[0]][j]]==(-1))
            {
               flags[edge[nn[0]][j]]=cluster;
               tmp++;
               nn[tmp]=edge[nn[0]][j];
               sizecl[cluster]++;
            }
         }

         nn[0]=nn[tmp];
         nn[tmp]=(-1);
         tmp--;
         if(nn[0]==(-1))   {break;}
      }

      r_node+=sizecl[cluster];
   }
   if(r_node==n_node)   {break;}

}

for(i=0,tmp=0;i<=cluster;i++)
{
   if(sizecl[i]>tmp)
   {
      largest_cluster=i;
      tmp=sizecl[i];
   }
}


free(nn);
free(sizecl);

return largest_cluster;

//return cluster;(Percolation);

}

void BA_network(int**edge, int*linkn, int n_node, int init_node, int add_link)	{

int i,j;
int current_node;
int dep,des;
int active_list[n_node];
int active,ch_active;
int add;
int add_node[add_link];

double prob,sum_linkn,prob_linkn;


for(i=0;i<n_node;i++)	{
	edge[i]=(int*)realloc(edge[i],sizeof(int)*1);
	edge[i][0]=(-1);
	linkn[i]=0;
}

current_node=init_node;

for(i=0;i<init_node-1;i++)	{
	for(j=i+1;j<init_node;j++)	{
		dep=i;
		des=j;

		linkn[dep]++;
		linkn[des]++;

		edge[dep]=(int*)realloc(edge[dep],sizeof(int)*linkn[dep]);
		edge[des]=(int*)realloc(edge[des],sizeof(int)*linkn[des]);

		edge[dep][linkn[dep]-1]=des;
		edge[des][linkn[des]-1]=dep;
	}
}

while(1)	{
	active=current_node;
	for(i=0;i<current_node;i++)	{
		active_list[i]=i;
	}
	add=0;
	while(1)	{
/*
		prob=ran4();
		for(i=0,sum_linkn=0.;i<active;i++)	{
			sum_linkn+=linkn[active_list[i]];
		}

		for(i=0,prob_linkn=0.;i<active;i++)	{
			prob_linkn+=(double)linkn[active_list[i]]/(double)sum_linkn;
			if(prob_linkn>prob)	{
				add_node[add]=active_list[i];
				break;
			}
		}
		active_list[i]=active_list[active-1];
*/
		ch_active=(int)(ran4()*active);
		add_node[add]=active_list[ch_active];
		active_list[ch_active]=active_list[active-1];

		active--;
		add++;
		if(add==add_link)	break;
	}

	for(i=0;i<add_link;i++)	{
		dep=current_node;
		des=add_node[i];

		linkn[dep]++;
		linkn[des]++;

		edge[dep]=(int*)realloc(edge[dep],sizeof(int)*linkn[dep]);
		edge[des]=(int*)realloc(edge[des],sizeof(int)*linkn[des]);

		edge[dep][linkn[dep]-1]=des;
		edge[des][linkn[des]-1]=dep;
	}

	current_node++;

	if(current_node==n_node)	break;
}

}

		



double sokolov_rewiring_method(int** net, int* linkn, int* flags, int largest_cluster, int n_node, double REWIRING_PROB, int CHECK)	{

int i,j;
int n_link;
int *linkarrayfrom,*linkarrayto;
int oldnode[4],oldnei[4],orderednode[4],active_list[4];
int link1,link2,tmp,active,ch_active,step;
double pearson_coefficient;


step=1;
n_link=0;

for(i=0;i<n_node;i++)	{
	for(j=0;j<linkn[i];j++)	{
		if(net[i][j]>i)   {n_link++;}
	}
}

linkarrayfrom=(int*)malloc(sizeof(int)*n_link);
linkarrayto=(int*)malloc(sizeof(int)*n_link);

n_link=0;

for(i=0;i<n_node;i++)	{
  for(j=0;j<linkn[i];j++)   {
		if(net[i][j]>i)	{
			linkarrayfrom[n_link]=i;
			linkarrayto[n_link]=net[i][j];
		
			n_link++;
		}
	}
}

while(1)	{
	link1=(int)(ran4()*n_link);
	while(1)	{
		link2=(int)(ran4()*n_link);
		if(link1!=link2)	{
			if(linkarrayfrom[link1]!=linkarrayfrom[link2] && linkarrayfrom[link1]!=linkarrayto[link2])	{
				if(linkarrayto[link1]!=linkarrayfrom[link2] && linkarrayto[link1]!=linkarrayto[link2]) {
					break;
				}
			}
		}
	}
	oldnode[0]=linkarrayfrom[link1];
	oldnode[1]=linkarrayto[link1];
	oldnode[2]=linkarrayfrom[link2];
	oldnode[3]=linkarrayto[link2];

	if(REWIRING_PROB>ran4())	{
		orderednode[0]=oldnode[0];

		if(linkn[oldnode[1]]>linkn[orderednode[0]])	{
			orderednode[1]=orderednode[0];
			orderednode[0]=oldnode[1];
		}
		else orderednode[1]=oldnode[1];

		if(linkn[oldnode[2]]>linkn[orderednode[0]])	{
			orderednode[2]=orderednode[1];
			orderednode[1]=orderednode[0];
			orderednode[0]=oldnode[2];
		}
		else	{
			if(linkn[oldnode[2]]>linkn[orderednode[1]])	{
				orderednode[2]=orderednode[1];
				orderednode[1]=oldnode[2];
			}
			else orderednode[2]=oldnode[2];
		}

		if(linkn[oldnode[3]]>linkn[orderednode[0]])	{
			orderednode[3]=orderednode[2];
			orderednode[2]=orderednode[1];
			orderednode[1]=orderednode[0];
			orderednode[0]=oldnode[3];
		}
		else	{
			if(linkn[oldnode[3]]>linkn[orderednode[1]])	{
				orderednode[3]=orderednode[2];
				orderednode[2]=orderednode[1];
				orderednode[1]=oldnode[3];
			}
			else	{
				if(linkn[oldnode[3]]>linkn[orderednode[2]])	{
					orderednode[3]=orderednode[2];
					orderednode[2]=oldnode[3];
				}
				else	{
					orderednode[3]=oldnode[3];
				}
			}
		}
	}
	
	else	{
		active=4;
		for(i=0;i<4;i++)  {active_list[i]=i;}
		for(i=0;i<4;i++)      {
			ch_active=ran4()*active;
			orderednode[i]=oldnode[ch_active];

			active_list[ch_active]=active_list[active-1];
			active--;
		}
	}



for(i=0,tmp=0;i<linkn[orderednode[1]];i++)	{
	if(net[orderednode[1]][i]==orderednode[0])   {tmp=1;break;}
}

if(tmp==0)	{
	for(i=0;i<linkn[orderednode[3]];i++)	{
		if(net[orderednode[3]][i]==orderednode[2])	{tmp=1;break;}
	}
}

if(tmp==0)	{
	step++;
	for(i=0;i<4;i++)	{
		for(j=0;j<4;j++)	{
			if(oldnode[j]==orderednode[i])	{
				switch(j)	{
					case 0 : oldnei[i]=oldnode[1];   break;
					case 1 : oldnei[i]=oldnode[0];   break;
					case 2 : oldnei[i]=oldnode[3];   break;
					case 3 : oldnei[i]=oldnode[2];   break;
				}
				break;
			}
		}
	}

	for(i=0;i<linkn[orderednode[0]];i++)	{
		if(net[orderednode[0]][i]==oldnei[0])	{
			net[orderednode[0]][i]=orderednode[1];
			break;
		}
	}
	for(i=0;i<linkn[orderednode[1]];i++)	{
		if(net[orderednode[1]][i]==oldnei[1])	{
			net[orderednode[1]][i]=orderednode[0];
			break;
		}
	}
	for(i=0;i<linkn[orderednode[2]];i++)	{
		if(net[orderednode[2]][i]==oldnei[2])	{
			net[orderednode[2]][i]=orderednode[3];
			break;
		}
	}
	for(i=0;i<linkn[orderednode[3]];i++)	{
		if(net[orderednode[3]][i]==oldnei[3])	{
			net[orderednode[3]][i]=orderednode[2];
			break;
		}
	}

	linkarrayfrom[link1]=orderednode[0];
	linkarrayto[link1]=orderednode[1];
	linkarrayfrom[link2]=orderednode[2];
	linkarrayto[link2]=orderednode[3];

}

if(step%CHECK==0)	{
	largest_cluster=Burning_algorithm(net,linkn,flags,n_node);
	pearson_coefficient=degree_correlation(net,linkn,flags,largest_cluster,n_node);
	printf("%f\n",pearson_coefficient);

	break;
}
}

for(i=0,tmp=0;i<n_node;i++)	{
	if(flags[i]==largest_cluster)	{tmp++;}
}

printf("real_node=%d\n",tmp);

free(linkarrayfrom);
free(linkarrayto);

return pearson_coefficient;

}

double degree_correlation(int** net, int* linkn, int* flags, int largest_cluster,int n_node)
{

int i,j;
double pearson_coefficient;
double x1,x2,x3;
double tot_link;

tot_link=0;
x1=x2=x3=0.;

for(i=0;i<n_node;i++)
{
   if(flags[i]==largest_cluster)
   {
      for(j=0;j<linkn[i];j++)
      {
         x1+=(double)linkn[i]*(double)linkn[net[i][j]];
         x2+=0.5*((double)linkn[i]+(double)linkn[net[i][j]]);
         x3+=0.5*((double)linkn[i]*(double)linkn[i]+(double)linkn[net[i][j]]*(double)linkn[net[i][j]]);
      }
      tot_link+=(double)linkn[i];
   }
}
tot_link/=2.0;

x1/=(double)(tot_link*2.0);
x2/=(double)(tot_link*2.0);   x2*=x2;
x3/=(double)(tot_link*2.0);

pearson_coefficient=(double)(x1-x2)/(double)(x3-x2);

return pearson_coefficient;

}


void pk_distribution(int *linkn, int *flags, int *pk, int n_node, int largest_cluster)
{

int i;

for(i=0;i<n_node;i++)	{pk[i]=0;}


for(i=0;i<n_node;i++)
{
   if(flags[i]==largest_cluster)
   {
      pk[linkn[i]]++;
   }
}

for(i=0;i<n_node;i++)
{
   if(pk[i]!=0)
   {
      printf("pk[%d]=%d\n",i,pk[i]);
   }
}


}

void knn_distribution(int **net, int *linkn, int *flags, double *knn, int *N_knn, int n_node)
{

int i,j,tmp;

for(i=0;i<n_node;i++)
{
	knn[i]=N_knn[i]=0;
}

for(i=0;i<n_node;i++)
{
   for(j=0,tmp=0;j<linkn[i];j++)
   {
      tmp+=linkn[net[i][j]];
   }

   knn[linkn[i]]+=(double)tmp/(double)linkn[i];
   N_knn[linkn[i]]++;

}

}

double MST_algorithm(int** edge, int* linkn, int* flags, int n_node, int** MST, int* MST_linkn)	{

int i,j,k;
int *active_list,*active_node;
int r_node=0;
int link_from,link_to;
int largest_cluster,active,act_node;
int tmp;
int dep,des;

double **weight;
double weight_value;
double min_weight,value;
double pearson_coefficient;

active_list=(int*)malloc(sizeof(int)*n_node);
active_node=(int*)malloc(sizeof(int)*n_node);

weight=(double**)malloc(sizeof(double*)*n_node);

for(i=0;i<n_node;i++)
{
	MST[i]=(int*)realloc(MST[i],sizeof(int));
	MST_linkn[i]=0;
	weight[i]=(double*)malloc(sizeof(double)*linkn[i]);
	MST[i][0]=active_list[i]=MST_linkn[i]=0;
	active_node[i]=(-1);
	for(j=0;j<linkn[i];j++)	weight[i][j]=0;
}

largest_cluster=Burning_algorithm(edge,linkn,flags,n_node);

for(i=0,r_node=0;i<n_node;i++)  {
	if(flags[i]==largest_cluster)   {
		r_node++;
		for(j=0;j<linkn[i];j++) {
			if(i>edge[i][j])    {
				weight_value=1./((double)linkn[i]*(double)linkn[edge[i][j]]);
//				weight_value=(double)linkn[i]*(double)linkn[edge[i][j]];
//				weight_value=ran4();
				weight[i][j]=weight_value;
				dep=edge[i][j];
 
				for(k=0;k<linkn[dep];k++)   {
					if(edge[dep][k]==i) {
						weight[dep][k]=weight_value;
						break;
					}
				}
			}
		}
	}
}
 
for(i=0,active=0;i<n_node;i++)  {
	if(flags[i]==largest_cluster)   {
		active=1;
		active_node[active-1]=i;
 
		break;
	}
}

for(i=0;i<n_node;i++)   {
	active_list[i]=0;
}
 
active_list[active_node[0]]=1;
 
while(1)    {
	min_weight=(double)INT_MAX;
	for(i=0;i<active;i++)   {
		dep=active_node[i];
		for(j=0;j<linkn[dep];j++)   {
			des=edge[dep][j];
			if(active_list[des]==0) {
				value=weight[dep][j];
				if(value<min_weight)    {
					link_from=dep;
					link_to=des;
					
					min_weight=value;
				}
			}
		}
	}

	active_list[link_to]=1;

	active++;
	active_node[active-1]=link_to;
//	if(active%1000==0)	printf("active=%d\n",active);

	MST_linkn[link_from]++;
	MST_linkn[link_to]++;

	MST[link_from]=(int*)realloc(MST[link_from],sizeof(int)*MST_linkn[link_from]);
	MST[link_to]=(int*)realloc(MST[link_to],sizeof(int)*MST_linkn[link_to]);

	MST[link_from][MST_linkn[link_from]-1]=link_to;
	MST[link_to][MST_linkn[link_to]-1]=link_from;

	if(r_node==active)  break;

}

pearson_coefficient=degree_correlation(MST,MST_linkn,flags,largest_cluster,n_node);

for(i=0;i<n_node;i++)	{
    free(weight[i]);
}

free(weight);
free(active_list);
free(active_node);


return pearson_coefficient;

}

double sokolov_rewiring_method_disasso(int** net, int* linkn, int* flags, int largest_cluster, int n_node, double REWIRING_PROB, int CHECK)	{

int i,j;
int n_link;
int *linkarrayfrom,*linkarrayto;
int oldnode[4],oldnei[4],orderednode[4],active_list[4];
int link1,link2,tmp,active,ch_active,step;
double pearson_coefficient;


step=1;
n_link=0;

for(i=0;i<n_node;i++)	{
	for(j=0;j<linkn[i];j++)	{
		if(net[i][j]>i)   {n_link++;}
	}
}

linkarrayfrom=(int*)malloc(sizeof(int)*n_link);
linkarrayto=(int*)malloc(sizeof(int)*n_link);

n_link=0;

for(i=0;i<n_node;i++)	{
  for(j=0;j<linkn[i];j++)   {
		if(net[i][j]>i)	{
			linkarrayfrom[n_link]=i;
			linkarrayto[n_link]=net[i][j];
		
			n_link++;
		}
	}
}

while(1)	{
	link1=(int)(ran4()*n_link);
	while(1)	{
		link2=(int)(ran4()*n_link);
		if(link1!=link2)	{
			if(linkarrayfrom[link1]!=linkarrayfrom[link2] && linkarrayfrom[link1]!=linkarrayto[link2])	{
				if(linkarrayto[link1]!=linkarrayfrom[link2] && linkarrayto[link1]!=linkarrayto[link2]) {
					break;
				}
			}
		}
	}
	oldnode[0]=linkarrayfrom[link1];
	oldnode[1]=linkarrayto[link1];
	oldnode[2]=linkarrayfrom[link2];
	oldnode[3]=linkarrayto[link2];

	if(REWIRING_PROB>ran4())	{
		orderednode[0]=oldnode[0];

		if(linkn[oldnode[1]]>linkn[orderednode[0]])	{
			orderednode[1]=orderednode[0];
			orderednode[0]=oldnode[1];
		}
		else orderednode[1]=oldnode[1];

		if(linkn[oldnode[2]]>linkn[orderednode[0]])	{
			orderednode[2]=orderednode[1];
			orderednode[1]=orderednode[0];
			orderednode[0]=oldnode[2];
		}
		else	{
			if(linkn[oldnode[2]]>linkn[orderednode[1]])	{
				orderednode[2]=orderednode[1];
				orderednode[1]=oldnode[2];
			}
			else orderednode[2]=oldnode[2];
		}

		if(linkn[oldnode[3]]>linkn[orderednode[0]])	{
			orderednode[3]=orderednode[2];
			orderednode[2]=orderednode[1];
			orderednode[1]=orderednode[0];
			orderednode[0]=oldnode[3];
		}
		else	{
			if(linkn[oldnode[3]]>linkn[orderednode[1]])	{
				orderednode[3]=orderednode[2];
				orderednode[2]=orderednode[1];
				orderednode[1]=oldnode[3];
			}
			else	{
				if(linkn[oldnode[3]]>linkn[orderednode[2]])	{
					orderednode[3]=orderednode[2];
					orderednode[2]=oldnode[3];
				}
				else	{
					orderednode[3]=oldnode[3];
				}
			}
		}
	}
	
	else	{
		active=4;
		for(i=0;i<4;i++)  {active_list[i]=i;}
		for(i=0;i<4;i++)      {
			ch_active=ran4()*active;
			orderednode[i]=oldnode[ch_active];

			active_list[ch_active]=active_list[active-1];
			active--;
		}
	}



for(i=0,tmp=0;i<linkn[orderednode[3]];i++)	{
	if(net[orderednode[3]][i]==orderednode[0])   {tmp=1;break;}
}

if(tmp==0)	{
	for(i=0;i<linkn[orderednode[1]];i++)	{
		if(net[orderednode[1]][i]==orderednode[2])	{tmp=1;break;}
	}
}

if(tmp==0)	{
	step++;
	for(i=0;i<4;i++)	{
		for(j=0;j<4;j++)	{
			if(oldnode[j]==orderednode[i])	{
				switch(j)	{
					case 0 : oldnei[i]=oldnode[1];   break;
					case 1 : oldnei[i]=oldnode[0];   break;
					case 2 : oldnei[i]=oldnode[3];   break;
					case 3 : oldnei[i]=oldnode[2];   break;
				}
				break;
			}
		}
	}

	for(i=0;i<linkn[orderednode[0]];i++)	{
		if(net[orderednode[0]][i]==oldnei[0])	{
			net[orderednode[0]][i]=orderednode[3];
			break;
		}
	}
	for(i=0;i<linkn[orderednode[1]];i++)	{
		if(net[orderednode[1]][i]==oldnei[1])	{
			net[orderednode[1]][i]=orderednode[2];
			break;
		}
	}
	for(i=0;i<linkn[orderednode[2]];i++)	{
		if(net[orderednode[2]][i]==oldnei[2])	{
			net[orderednode[2]][i]=orderednode[1];
			break;
		}
	}
	for(i=0;i<linkn[orderednode[3]];i++)	{
		if(net[orderednode[3]][i]==oldnei[3])	{
			net[orderednode[3]][i]=orderednode[0];
			break;
		}
	}

	linkarrayfrom[link1]=orderednode[0];
	linkarrayto[link1]=orderednode[3];
	linkarrayfrom[link2]=orderednode[1];
	linkarrayto[link2]=orderednode[2];

}

if(step%CHECK==0)	{
	largest_cluster=Burning_algorithm(net,linkn,flags,n_node);
	pearson_coefficient=degree_correlation(net,linkn,flags,largest_cluster,n_node);
	printf("%f\n",pearson_coefficient);

	break;
}
}

for(i=0,tmp=0;i<n_node;i++)	{
	if(flags[i]==largest_cluster)	{tmp++;}
}

printf("real_node=%d\n",tmp);

free(linkarrayfrom);
free(linkarrayto);

return pearson_coefficient;

}

double degree_correlation_modify(int** net, int* linkn, int n_node)
{

int i,j;
double pearson_coefficient;
double x1,x2,x3;
double tot_link;

tot_link=0;
x1=x2=x3=0.;

for(i=0;i<n_node;i++)
{
   if(linkn[i]!=0)
   {
      for(j=0;j<linkn[i];j++)
      {
         x1+=(double)linkn[i]*(double)linkn[net[i][j]];
         x2+=0.5*((double)linkn[i]+(double)linkn[net[i][j]]);
         x3+=0.5*((double)linkn[i]*(double)linkn[i]+(double)linkn[net[i][j]]*(double)linkn[net[i][j]]);
      }
      tot_link+=(double)linkn[i];
   }
}
tot_link/=2.0;

x1/=(double)(tot_link*2.0);
x2/=(double)(tot_link*2.0);   x2*=x2;
x3/=(double)(tot_link*2.0);

pearson_coefficient=(double)(x1-x2)/(double)(x3-x2);

return pearson_coefficient;

}
