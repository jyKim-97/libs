void SF_configuration(int*, int, double, double);

void SF_static(int**, int*, int, double, double);
//(network array, link number, node number, average_link, gamma)

void Random_ErdosRenyi(int**, int*, int, double);
//(network array, link number, node number, average_link)

void BA_network(int**, int*, int, int, int);

int Burning_algorithm(int**, int*, int*, int);
//(network array, link number, flags, node number)

int Burning_algorithm_modify(int**, int*, int*, int*, int);
//(network array, link number, flags, cluster_size,node number)

void n_dimension_lattice(int**, int*, int, int, int);
//(network array, link number, node number, length, dimension)

void n_dimension_tree(int**, int*, int, int, int);
//(network array, link number, node number, length, dimension)

void n_dimension_lattice_nb(int**, int*, int, int, int);
//(network array, link number, node number, length, dimension)

void pk_distribution(int*, int*, int *, int, int);
//(link number, flags, pk array, node number, cluster_size)

void knn_distribution(int**, int*, int*, double*, int*, int);
//(network array, link number, knn array, N_knn array, node number)

double sokolov_rewiring_method(int**, int*, int*, int, int, double, int);
//(network array, link number, flags, largest_cluster, node number)

double sokolov_rewiring_method_asso(int**, int*, int*, int, int, double, int);
//(network array, link number, flags, largest_cluster, node number)

double sokolov_rewiring_method_disasso(int**, int*, int*, int, int, double, int);
//(network array, link number, flags, largest_cluster, node number)

double degree_correlation(int**, int*, int*, int, int);
//(network array, link number, flags, largest_cluster, node number)

int Hoshen_Kopelman(int**, int*, int*, int);

int Hoshen_Kopelman_modify(int**, int*, int*, int, int, int);

double degree_correlation_modify(int**, int*, int);
//double MST_algorithm(int**, int*, int*, int, int**, int*);

//int clustercoloring2_forstatic(int*, int**, int*, int);
