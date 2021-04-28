#include <stdio.h>
#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

int n_knots = 0, n_elem = 0, n_cond = 0;

struct GlobalMatrix
{
	vector<double> ggl, ggu, di;
	vector<double> l, u, d;
	vector<double> r, p, z;
	vector<double> x, temp, temp0, F;
	vector<int> ig, jg;

};

struct knot
{
	double x;
	double y;
};

struct elem
{
	vector<knot> knts;
	vector<int> global_num;
	int lambda, func;
	double gamma;
};

struct localMatrix
{
	int elemNum;
	vector<vector<double>> localM;
};

struct condition
{
	int type;
	int vertex;
	int num;
};

double Lambda(int num, knot& k)
{
	switch (num)
	{
	case 1:
		return 1;
	case 2:
		return 1;
	default:
		break;
	}
	return 1;
}

double Func(int num, knot& k, elem &e)
{
	switch (num)
	{
	case(1):
		return k.x * k.x - 2;
	case(2):
		return k.x * k.x * k.x - 6 * k.x;
	case 3:
		return k.x * k.x * k.x * k.x - 12 * k.x * k.x;
	case 4:
		return sin(k.x * k.y) * (k.x * k.x + k.y * k.y + 1) + k.x * k.x * k.y - 2 * k.y;
	case 5:
		return 1;
	case 6:
		return k.x;
	case 7:
		return k.x * k.x * k.y - 2 * k.y;
	default:
		break;
	}
	
}

double FirstCondition(int num, knot& k)
{
	switch (num)
	{
	case 1:
		return k.x * k.x;
	case 2:
		return k.x * k.x * k.x;
	case 3:
		return k.x * k.x * k.x * k.x;
	case 4:
		return k.x * k.x * k.y + sin(k.x * k.y);
	case 5:
		return 1;
	case 6:
		return k.x;
	case 7:
		return k.x * k.x * k.y;
	default:
		break;
	}
}

bool Input(vector<knot>& knots, vector<elem> &elems, vector<condition> &conds)
{
	FILE* info;
	if (!fopen_s(&info, "info.txt", "r"))
	{
		fscanf_s(info, "%d %d", &n_knots, &n_elem);
	}
	else
	{
		cout << "info.txt is invalid";
		return false;
	}
	fclose(info);

	knots.resize(n_knots);
	elems.resize(n_elem);

	FILE* xy;
	if (!fopen_s(&xy, "xy.txt", "r"))
	{
		for (int i = 0; i < n_knots; i++)
		{
			fscanf_s(xy, "%lf %lf", &(knots[i].x), &(knots[i].y));
		}
	}
	else
	{
		cout << "xy.txt is invalid";
		return false;
	}
	fclose(xy);

	FILE* el;
	if (!fopen_s(&el, "elem.txt", "r"))
	{
		for (int i = 0; i < n_elem; i++)
		{
			elems[i].global_num.resize(6);
			elems[i].knts.resize(6);
			fscanf_s(el, "%d %d %d %d %d %d ", &(elems[i].global_num[0]), &(elems[i].global_num[1]),
				&(elems[i].global_num[2]), &(elems[i].global_num[3]), &(elems[i].global_num[4]), &(elems[i].global_num[5]));
			for (int j = 0; j < 6; j++)
			{
				elems[i].global_num[j]--;
				elems[i].knts[j] = knots[elems[i].global_num[j]];
			}
		}
	}
	else
	{
		cout << "elem.txt is invalid";
		return false;
	}
	fclose(el);

	FILE* mat;
	if (!fopen_s(&mat, "mat.txt", "r"))
	{
		for (int i = 0; i < n_elem; i++)
		{
			fscanf_s(mat, "%d %d %lf", &(elems[i].lambda), &(elems[i].func), &(elems[i].gamma));
		}
	}
	else
	{
		cout << "mat.txt is invalid";
		return false;
	}
	fclose(mat);

	FILE* s;
	if (!fopen_s(&s, "s1.txt", "r"))
	{
		fscanf_s(s, "%d", &n_cond);
		conds.resize(n_cond);
		for (int i = 0; i < n_cond; i++)
		{
			fscanf_s(s, "%d %d %d", &(conds[i].type), &(conds[i].num), &(conds[i].vertex));
		}
	}
	else
	{
		cout << "s1.txt is invalid";
		return false;
	}
	fclose(s);

	return true;
}

double DetD(elem &el)
{
	double det = (el.knts[1].x - el.knts[0].x) * (el.knts[2].y - el.knts[0].y) -
		(el.knts[2].x - el.knts[0].x) * (el.knts[1].y - el.knts[0].y);
	return det;
}

void CreateAlphas(vector<vector<double>> &alphas, elem &el, double &detD)
{
	alphas[0][0] = (el.knts[1].y - el.knts[2].y) / detD;
	alphas[0][1] = (el.knts[2].x - el.knts[1].x) / detD;
	alphas[1][0] = (el.knts[2].y - el.knts[0].y) / detD;
	alphas[1][1] = (el.knts[0].x - el.knts[2].x) / detD;
	alphas[2][0] = (el.knts[0].y - el.knts[1].y) / detD;
	alphas[2][1] = (el.knts[1].x - el.knts[0].x) / detD;
}

void addElemToGlobal(GlobalMatrix& M, int& i, int& j, double& value)
{
	int ind;

	if (i == j)
	{
		M.di[i] += value;
	}
	else
	{
		if (i > j)
		{
			for (ind = M.ig[i]; ind < M.ig[i + 1]; ind++)
			{
				if (M.jg[ind] == j)
				{
					M.ggl[ind] += value;
					break;
				}
			}
		}
		else
		{
			for (ind = M.ig[j]; ind < M.ig[j + 1]; ind++)
			{
				if (M.jg[ind] == i)
				{
					M.ggu[ind] += value;
					break;
				}
			}
		}
	}
}

void addLocalToGlobal(GlobalMatrix& M, elem& el, vector<vector<double>>& localM, vector<double> &localB)
{
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			addElemToGlobal(M, el.global_num[i], el.global_num[j], localM[i][j]);
		}
	}

	for (int i = 0; i < 6; i++)
	{
		M.F[el.global_num[i]] = localB[i];
	}
}

void CreateArrayOfLocals(vector<localMatrix> &arrayLocals, vector<elem> &elems)
{
	vector<vector<double>> M(6), G(6), alphas(3);
	double detD = 0, buff = 0., l1 = 0., l2 = 0., l3 = 0.;
	for (int i = 0; i < 6; i++)
	{
		M[i].resize(6);
		G[i].resize(6);
	}

	for (int i = 0; i < 3; i++)
	{
		alphas[i].resize(2);
	}

	for (int i = 0; i < n_elem; i++)
	{
		arrayLocals[i].elemNum = i;
		arrayLocals[i].localM.resize(6);
		for (int j = 0; j < 6; j++)
		{
			arrayLocals[i].localM[j].resize(6);
		}

		detD = DetD(elems[i]);
		CreateAlphas(alphas, elems[i], detD);

		detD = abs(detD);
		buff = elems[i].gamma;
		M[0][0] = (1. / 60) * buff;
		M[1][1] = (1. / 60) * buff;
		M[2][2] = (1. / 60) * buff;
		M[3][3] = (4. / 45) * buff;
		M[4][4] = (4. / 45) * buff;
		M[5][5] = (4. / 45) * buff;
		M[0][1] = M[1][0] = (-1. / 360) * buff;
		M[0][2] = M[2][0] = (-1. / 360) * buff;
		M[1][2] = M[2][1] = (-1. / 360) * buff;
		M[3][4] = M[4][3] = (2. / 45) * buff;
		M[3][5] = M[5][3] = (2. / 45) * buff;
		M[4][5] = M[5][4] = (2. / 45) * buff;
		M[2][3] = M[3][2] = M[0][4] = M[4][0] = M[1][5] = M[5][1] = (-1. / 90) * buff;

		l1 = Lambda(elems[i].lambda, elems[i].knts[0]);
		l2 = Lambda(elems[i].lambda, elems[i].knts[1]);
		l3 = Lambda(elems[i].lambda, elems[i].knts[2]);

		G[0][0] = (alphas[0][0] * alphas[0][0] + alphas[0][1] * alphas[0][1]) * (3 * l1 + l2 + l3) / 10;
		G[1][1] = (alphas[1][0] * alphas[1][0] + alphas[1][1] * alphas[1][1]) * (l1 + 3 * l2 + l3) / 10;
		G[2][2] = (alphas[2][0] * alphas[2][0] + alphas[2][1] * alphas[2][1]) * (l1 + l2 + 3 * l3) / 10;

		G[3][3] = (4. / 15) * ((alphas[0][0] * alphas[0][0] + alphas[0][1] * alphas[0][1]) * (l1 + 3 * l2 + l3) +
			(alphas[0][0] * alphas[1][0] + alphas[0][1] * alphas[1][1]) * (2 * l1 + 2 * l2 + l3) +
			(alphas[1][0] * alphas[1][0] + alphas[1][1] * alphas[1][1]) * (3 * l1 + l2 + l3));

		G[4][4] = (4. / 15) * ((alphas[1][0] * alphas[1][0] + alphas[1][1] * alphas[1][1]) * (l1 + l2 + 3 * l3) +
			(alphas[1][0] * alphas[2][0] + alphas[1][1] * alphas[2][1]) * (l1 + 2 * l2 + 2 * l3) +
			(alphas[2][0] * alphas[2][0] + alphas[2][1] * alphas[2][1]) * (l1 + 3 * l2 + l3));

		G[5][5] = (4. / 15) * ((alphas[0][0] * alphas[0][0] + alphas[0][1] * alphas[0][1]) * (l1 + l2 + 3 * l3) +
			(alphas[0][0] * alphas[2][0] + alphas[0][1] * alphas[2][1]) * (2 * l1 + l2 + 2 * l3) +
			(alphas[2][0] * alphas[2][0] + alphas[2][1] * alphas[2][1]) * (3 * l1 + l2 + l3));

		G[0][1] = G[1][0] = -(1. / 30) * (alphas[0][0] * alphas[1][0] + alphas[0][1] * alphas[1][1]) * (2 * l1 + 2 * l2 + l3);
		G[0][2] = G[2][0] = -(1. / 30) * (alphas[0][0] * alphas[2][0] + alphas[0][1] * alphas[2][1]) * (2 * l1 + l2 + 2 * l3);
		G[1][2] = G[2][1] = -(1. / 30) * (alphas[1][0] * alphas[2][0] + alphas[1][1] * alphas[2][1]) * (l1 + 2 * l2 + 2 * l3);

		G[3][4] = G[4][3] = (4. / 15) * ((alphas[0][0] * alphas[1][0] + alphas[0][1] * alphas[1][1]) * (l1 / 2 + l2 + l3) +
			(alphas[0][0] * alphas[2][0] + alphas[0][1] * alphas[2][1]) * (l1 + 3 * l2 + l3) +
			(alphas[1][0] * alphas[1][0] + alphas[1][1] * alphas[1][1]) * (l1 + l2 / 2 + l3) +
			(alphas[1][0] * alphas[2][0] + alphas[1][1] * alphas[2][1]) * (l1 + l2 + l3 / 2));

		G[3][5] = G[5][3] = (4. / 15) * ((alphas[0][0] * alphas[2][0] + alphas[0][1] * alphas[2][1]) * (l1 + l2 + l3 / 2) +
			(alphas[0][0] * alphas[1][0] + alphas[0][1] * alphas[1][1]) * (l1 + l2 / 2 + l3) +
			(alphas[0][0] * alphas[0][0] + alphas[0][1] * alphas[0][1]) * (l1 / 2 + l2 + l3) +
			(alphas[1][0] * alphas[2][0] + alphas[1][1] * alphas[2][1]) * (3 * l1 + l2 + l3));

		G[4][5] = G[5][4] = (4. / 15) * ((alphas[0][0] * alphas[2][0] + alphas[0][1] * alphas[2][1]) * (l1 / 2 + l2 + l3) +
			(alphas[1][0] * alphas[2][0] + alphas[1][1] * alphas[2][1]) * (l1 + l2 / 2 + l3) +
			(alphas[2][0] * alphas[2][0] + alphas[2][1] * alphas[2][1]) * (l1 + l2 + l3 / 2) +
			(alphas[0][0] * alphas[1][0] + alphas[0][1] * alphas[1][1]) * (l1 + l2 + 3 * l3));

		G[0][3] = G[3][0] = (1. / 30) * ((alphas[0][0] * alphas[0][0] + alphas[0][1] * alphas[0][1]) * (3 * l1 - 2 * l2 - l3) +
			(alphas[0][0] * alphas[1][0] + alphas[0][1] * alphas[1][1]) * (14 * l1 + 3 * l2 + 3 * l3));

		G[1][3] = G[3][1] = (1. / 30) * ((alphas[1][0] * alphas[1][0] + alphas[1][1] * alphas[1][1]) * (-2 * l1 + 3 * l2 - l3) +
			(alphas[0][0] * alphas[1][0] + alphas[0][1] * alphas[1][1]) * (3 * l1 + 14 * l2 + 3 * l3));

		G[2][3] = G[3][2] = (1. / 30) * ((alphas[1][0] * alphas[2][0] + alphas[1][1] * alphas[2][1]) * (-2 * l1 - l2 + 3 * l3) +
			(alphas[0][0] * alphas[2][0] + alphas[0][1] * alphas[2][1]) * (-1 * l1 - 2 * l2 + 3 * l3));

		G[0][4] = G[4][0] = (1. / 30) * ((alphas[1][0] * alphas[0][0] + alphas[0][1] * alphas[1][1]) * (3 * l1 - l2 - 2 * l3) +
			(alphas[0][0] * alphas[2][0] + alphas[0][1] * alphas[2][1]) * (3 * l1 - 2 * l2 - l3));

		G[1][4] = G[4][1] = (1. / 30) * ((alphas[1][0] * alphas[0][0] + alphas[1][1] * alphas[0][1]) * (-1 * l1 + 3 * l2 - 2 * l3) +
			(alphas[1][0] * alphas[2][0] + alphas[1][1] * alphas[2][1]) * (3 * l1 + 14 * l2 + 3 * l3));

		G[2][4] = G[4][2] = (1. / 30) * ((alphas[1][0] * alphas[2][0] + alphas[1][1] * alphas[2][1]) * (3 * l1 + 3 * l2 + 14 * l3) +
			(alphas[2][0] * alphas[2][0] + alphas[2][1] * alphas[2][1]) * (-1 * l1 - 2 * l2 + 3 * l3));

		G[0][5] = G[5][0] = (1. / 30) * ((alphas[0][0] * alphas[0][0] + alphas[0][1] * alphas[0][1]) * (3 * l1 - l2 - 2 * l3) +
			(alphas[0][0] * alphas[2][0] + alphas[0][1] * alphas[2][1]) * (14 * l1 + 3 * l2 + 3 * l3));

		G[1][5] = G[5][1] = (1. / 30) * ((alphas[0][0] * alphas[1][0] + alphas[0][1] * alphas[1][1]) * (-1 * l1 + 3 * l2 - 2 * l3) +
			(alphas[1][0] * alphas[2][0] + alphas[1][1] * alphas[2][1]) * (-2 * l1 + 3 * l2 - l3));

		G[2][5] = G[5][2] = (1. / 30) * ((alphas[0][0] * alphas[2][0] + alphas[0][1] * alphas[2][1]) * (3 * l1 + 3 * l2 + 14 * l3) +
			(alphas[2][0] * alphas[2][0] + alphas[2][1] * alphas[2][1]) * (-2 * l1 - l2 + 3 * l3));

		for (int k = 0; k < 6; k++)
		{
			for (int z = 0; z < 6; z++)
			{
				arrayLocals[i].localM[k][z] = detD * (M[k][z] + G[k][z]);
			}
		}
	}
}

void CreateArrayOfLocalB(vector<vector<double>> &arr, vector<elem> &elems)
{
	double buff = 0., f1 = 0., f2 = 0., f3 = 0., f4 = 0., f5 = 0., f6 = 0.;
	for (int i = 0; i < n_elem; i++)
	{
		arr[i].resize(6);
		double detD = abs(DetD(elems[i]));
		buff = detD / 360;

		f1 = Func(elems[i].func, elems[i].knts[0], elems[i]);
		f2 = Func(elems[i].func, elems[i].knts[1], elems[i]);
		f3 = Func(elems[i].func, elems[i].knts[2], elems[i]);
		f4 = Func(elems[i].func, elems[i].knts[3], elems[i]);
		f5 = Func(elems[i].func, elems[i].knts[4], elems[i]);
		f6 = Func(elems[i].func, elems[i].knts[5], elems[i]);

		arr[i][0] = buff * (6 * f1 - f2 - f3 - 4 * f5);
		arr[i][1] = buff * (-f1 + f2 * 6 - f3 - 4 * f6);
		arr[i][2] = buff * (-f1 - f2 + 6 * f3 - 4 * f4);

		buff = detD / 45;

		arr[i][3] = buff * (-f3 / 2 + 4 * f4 + 2 * f5 + 2 * f6);
		arr[i][4] = buff * (-f1 / 2 + 2 * f4 + 4 * f5 + 2 * f6);
		arr[i][5] = buff * (-f2 / 2 + 2 * f4 + 2 * f5 + 4 * f6);
	}
}

void CreateGlobalSystem(vector<vector<double>> &A, vector<double> &b, vector<localMatrix>& arrA, vector<vector<double>>& arrB, vector<elem>& elems)
{
	for (int k = 0; k < n_elem; k++)
	{
		for (int i = 0; i < 6; i++)
		{
			int indexB = elems[k].global_num[i];
			int indexI = elems[arrA[k].elemNum].global_num[i];
			for (int j = 0; j < 6; j++)
			{
				int indexJ = elems[arrA[k].elemNum].global_num[j];
				A[indexI][indexJ] += arrA[k].localM[i][j];
			}
			b[indexB] += arrB[k][i];
		}
	}
}

void AccountConditions_1(vector<vector<double>> &A, vector<double> &b, vector<condition> &conds, vector<knot> &knots)
{
	for (int i = 0; i < n_cond; i++)
	{
		int index = conds[i].vertex - 1;
		for (int j = 0; j < n_knots; j++)
		{
			A[index][j] = (index == j ? 1 : 0);
		}
		b[index] = FirstCondition(conds[i].num, knots[index]);
	}
}

void AccountConditions(GlobalMatrix& A, vector<double>& b, vector<condition>& conds, vector<knot>& knots)
{
	int ibeg, iend, jind;
	for (int i = 0; i < n_cond; i++)
	{
		int index = conds[i].vertex - 1;
		A.di[index] = 1;
		b[index] = FirstCondition(conds[i].num, knots[index]);

		ibeg = A.ig[index];
		iend = A.ig[index + 1];
		for (int j = ibeg; j < iend; j++)
		{
			A.ggl[j] = 0;
		}
		// Обнуляем внедиагональные элементы i-ой строки в ggu
		jind = index;
		for (int j = 0; j < A.ig[n_knots]; j++)
		{
			if (A.jg[j] == jind)
			{
				A.ggu[j] = 0;
			}
		}
	}
}

double ScalarMult(vector<double> &vec1, vector<double>& vec2)
{
	double res = 0.;
	for (int i = 0; i < n_knots; i++)
	{
		res += vec1[i] * vec2[i];
	}
	return res;
}

void MatrixVectorMult(vector<vector<double>>& A, vector<double> &vec, vector<double> &res)
{
	for (int i = 0; i < n_knots; i++)
	{
		res[i] = ScalarMult(A[i], vec);
	}
}

void LOC(vector<vector<double>> &A, vector<double> &b, vector<double> &x)
{
	double residual = 0., residual1 = 0., alpha = 0., beta = 0., scMult = 0., eps = 1e-15;
	vector<double> r(n_knots), z(n_knots), p(n_knots);

	for (int i = 0; i < n_knots; i++)
	{
		b[i] -= ScalarMult(A[i], x);
	}
	r = z = b;
	MatrixVectorMult(A, z, p);
	residual = ScalarMult(r, r);

	for (int i = 1; i < 100000 && residual > eps && residual != residual1; i++)
	{
		residual1 = residual;
		scMult = ScalarMult(p, p);
		alpha = ScalarMult(p, r) / scMult;

		for (int j = 0; j < n_knots; j++)
		{
			x[j] += alpha * z[j];
			r[j] -= alpha * p[j];
		}

		MatrixVectorMult(A, r, b);

		beta = -ScalarMult(p, b) / scMult;

		for (int j = 0; j < n_knots; j++)
		{
			z[j] = r[j] + beta * z[j];
			p[j] = b[j] + beta * p[j];
		}

		residual -= alpha * alpha * scMult;
	}
}

// Процедура построения портрета матрицы
void Portrait(GlobalMatrix &M, vector<elem> &elems)
{
	vector<int> listbeg(n_knots);
	vector<int> list1(8 * n_knots);
	vector<int> list2(8 * n_knots);
	
	int ielem, i, k, j;
	int ind1, ind2, iaddr;
	int listsize = -1;

	for (ielem = 0; ielem < n_elem; ielem++)
	{
		for (i = 0; i < 6; i++)
		{
			k = elems[ielem].global_num[i];
			for (j = i + 1; j < 6; j++)
			{
				ind1 = k;
				ind2 = elems[ielem].global_num[j];
				if (ind2 < ind1)
				{
					ind1 = ind2;
					ind2 = k;
				}
				iaddr = listbeg[ind2];
				if (iaddr == 0)
				{
					listsize++;
					listbeg[ind2] = listsize;
					list1[listsize] = ind1;
					list2[listsize] = 0;
				}
				else
				{
					while (list1[iaddr] < ind1 && list2[iaddr] > 0)
					{
						iaddr = list2[iaddr];
					}
					if (list1[iaddr] > ind1)
					{
						listsize++;
						list1[listsize] = list1[iaddr];
						list2[listsize] = list2[iaddr];
						list1[iaddr] = ind1;
						list2[iaddr] = listsize;
					}
					else
					{
						if (list1[iaddr] < ind1)
						{
							listsize++;
							list2[iaddr] = listsize;
							list1[listsize] = ind1;
							list2[listsize] = 0;
						}
					}
				}
			}
		}
	}

	M.ig.resize(n_knots + 1);
	M.di.resize(n_knots);
	M.ggl.resize(list1.size());
	M.ggu.resize(list1.size());
	M.u.resize(list1.size());
	M.l.resize(list1.size());
	M.d.resize(n_knots);
	M.jg.resize(list1.size());
	M.x.resize(n_knots);
	M.temp.resize(n_knots);
	M.temp0.resize(n_knots);
	M.F.resize(n_knots);
	M.r.resize(n_knots);
	M.p.resize(n_knots);
	M.z.resize(n_knots);

	M.ig[0] = 0;

	for (i = 0; i < n_knots; i++)
	{
		M.ig[i + 1] = M.ig[i];
		iaddr = listbeg[i];
		while (iaddr != 0)
		{
			M.jg[M.ig[i + 1] + 1] = list1[iaddr];
			M.ig[i + 1]++;
			iaddr = list2[iaddr];
		}
	}
	for (i = 3; i < n_knots + 1; i++)
	{
		M.ig[i]++;
	}
}

// LU Факторизация
void CalcLU(GlobalMatrix& matrix) {
	double sumU, sumL, sumD;
	int n = n_knots;
	matrix.l = matrix.ggl;
	matrix.u = matrix.ggu;
	matrix.d = matrix.di;


	for (int i = 0; i < n; i++) {

		sumD = 0;

		int begI = matrix.ig[i];
		int endI = matrix.ig[i + 1];
		for (int igi = begI; igi < endI; igi++) {

			sumU = 0;
			sumL = 0;

			int Jindex = matrix.jg[igi];

			for (int igj = begI; igj < igi; igj++) {

				int begJ = matrix.ig[Jindex];
				int endJ = matrix.ig[Jindex + 1];

				for (int jgi = begJ; jgi < endJ; jgi++) {

					if (matrix.jg[igj] == matrix.jg[jgi]) {

						sumL += matrix.l[igj] * matrix.u[jgi];

						sumU += matrix.l[jgi] * matrix.u[igj];

					}

				}
			}
			matrix.l[igi] -= sumL;
			matrix.u[igi] -= sumU;
			matrix.u[igi] /= matrix.d[Jindex];
			sumD += matrix.l[igi] * matrix.u[igi];
		}

		matrix.d[i] -= sumD;
	}
}

// Прямой ход Ly = F
void CalcDir(GlobalMatrix& matrix, vector<double>& y, vector<double>& F) {
	double sum, buf;
	int n = n_knots;


	for (int i = 0; i < n; i++) {
		y[i] = F[i];
	}

	for (int i = 0; i < n; i++) {

		sum = 0;

		int begI = matrix.ig[i];
		int endI = matrix.ig[i + 1];

		for (int igi = begI; igi < endI; igi++) {

			sum += y[matrix.jg[igi]] * matrix.l[igi];

		}

		buf = y[i] - sum;

		y[i] = buf / matrix.d[i];

	}

}

// Обратный ход Ux = y
void CalcRev(GlobalMatrix& matrix, vector<double>& x, vector<double>& y) {
	int n = n_knots;

	for (int i = 0; i < n; i++) {
		x[i] = y[i];
	}

	for (int i = n_knots - 1; i >= 0; i--) {

		int begI = matrix.ig[i];
		int endI = matrix.ig[i + 1];

		for (int igi = begI; igi < endI; igi++) {

			x[matrix.jg[igi]] -= x[i] * matrix.u[igi];

		}

	}

}

// Умножение матрицы на вектор Ax = res
void MultMV(GlobalMatrix& matrix, vector<double>& x, vector<double>& res) {
	int n = n_knots;

	for (int i = 0; i < n; i++) {

		res[i] = matrix.di[i] * x[i];

		int begI = matrix.ig[i];
		int endI = matrix.ig[i + 1];

		for (int igi = begI; igi < endI; igi++) {

			int Jindex = matrix.jg[igi];

			res[i] += matrix.ggl[igi] * x[Jindex];
			res[Jindex] += matrix.ggu[igi] * x[i];

		}
	}
}

// Скалярное произведение двух векторов
double ScalarProd(vector<double>& x, vector<double>& y) {
	int n = x.size();

	double result = 0;

	for (int i = 0; i < n; i++) {
		result += x[i] * y[i];
	}

	return result;
}

// Локально-оптимальная схема c факторизацией LU
void LOS_LU(GlobalMatrix& matrix) {

	double alpha, beta, norm, temp_nev = 0;

	int n = n_knots, maxiter = 1e9;
	double epsilon = 1e-15;

	CalcLU(matrix);
	// A * x0
	MultMV(matrix, matrix.x, matrix.temp);

	// f - A * x0
	for (int i = 0; i < n; i++) {
		matrix.temp[i] = matrix.F[i] - matrix.temp[i];
	}

	// L * r0 = f - A * x0
	CalcDir(matrix, matrix.r, matrix.temp);

	// U * z0 = r0
	CalcRev(matrix, matrix.z, matrix.r);

	// A * z0
	MultMV(matrix, matrix.z, matrix.temp);

	// L * p0 = A * z0
	CalcDir(matrix, matrix.p, matrix.temp);

	norm = ScalarProd(matrix.r, matrix.r);

	int k;

	for (k = 0; k < maxiter && norm > epsilon && norm != temp_nev; k++) {

		// если невязка не изменилась, то выходим из итерационного процесса
		temp_nev = norm;

		alpha = ScalarProd(matrix.p, matrix.r) / ScalarProd(matrix.p, matrix.p);

		for (int i = 0; i < n; i++) {
			matrix.x[i] = matrix.x[i] + alpha * matrix.z[i];
			matrix.r[i] = matrix.r[i] - alpha * matrix.p[i];
		}

		// U * temp = r
		CalcRev(matrix, matrix.temp, matrix.r);

		// A * U-1 * r = temp0
		MultMV(matrix, matrix.temp, matrix.temp0);

		// L * temp = A * U-1 * r 
		CalcDir(matrix, matrix.temp, matrix.temp0);

		beta = -1 * ScalarProd(matrix.p, matrix.temp) / ScalarProd(matrix.p, matrix.p);

		// U * temp0 = r
		CalcRev(matrix, matrix.temp0, matrix.r);

		norm = norm - alpha * alpha * ScalarProd(matrix.p, matrix.p);

		for (int i = 0; i < n; i++) {

			matrix.z[i] = matrix.temp0[i] + beta * matrix.z[i];
			matrix.p[i] = matrix.temp[i] + beta * matrix.p[i];

		}

	}

}

void rndPoint(knot &k, elem &el, vector<double> &q)
{
	vector<vector<double>> alpha(3);
	for (int i = 0; i < 3; i++)
		alpha[i].resize(3);
	vector<double> L(3);
	vector<double> basis(6);

	double detD = DetD(el);
	double S23 = abs((el.knts[2].x - el.knts[1].x) * (k.y - el.knts[1].y) - (k.x - el.knts[1].x) * (el.knts[2].y - el.knts[1].y));
	double S31 = abs((el.knts[0].x - el.knts[2].x) * (k.y - el.knts[2].y) - (k.x - el.knts[2].x) * (el.knts[0].y - el.knts[2].y));
	double S12 = abs((el.knts[1].x - el.knts[0].x) * (k.y - el.knts[0].y) - (k.x - el.knts[0].x) * (el.knts[1].y - el.knts[0].y));

	L[0] = S23 / abs(detD);
	L[1] = S31 / abs(detD);
	L[2] = S12 / abs(detD);

	for (int i = 0; i < 3; i++)
		basis[i] = L[i] * (2. * L[i] - 1);
	basis[3] = 4. * L[0] * L[1];
	basis[4] = 4. * L[1] * L[2];
	basis[5] = 4. * L[0] * L[2];
	double sum = 0;
	// для вершин
	sum += q[el.global_num[0]] * basis[0];
	sum += q[el.global_num[1]] * basis[1];
	sum += q[el.global_num[2]] * basis[2];

	// для середин
	sum += q[el.global_num[3]] * basis[3];
	sum += q[el.global_num[4]] * basis[4];
	sum += q[el.global_num[5]] * basis[5];
	double resh = FirstCondition(3, k);
	cout << setprecision(3) << sum << "\t" << "\t" << resh << "\t" << "\t" << fabs(resh - sum) << "\n\n";
}

int main()
{
	GlobalMatrix M;
	vector<knot> knots;
	vector<elem> elems;
	vector<condition> conds;
	if (!Input(knots, elems, conds)) return 1;

	Portrait(M, elems);
	GlobalMatrix G = M, A = M;

	vector<localMatrix> arrayLocals(n_elem);
	CreateArrayOfLocals(arrayLocals, elems);

	vector<vector<double>> arrayLocalsB(n_elem);
	CreateArrayOfLocalB(arrayLocalsB, elems);
 
	vector<vector<double>> A_1(n_knots);
	vector<double> b(n_knots);
	for (int i = 0; i < n_knots; i++)
		A_1[i].resize(n_knots);

	//CreateGlobalSystem(A_1, b, arrayLocals, arrayLocalsB, elems);

	for (int i = 0; i < n_elem; i++)
	{
		addLocalToGlobal(A, elems[i], arrayLocals[i].localM, arrayLocalsB[i]);
	}

	AccountConditions(A, A.F, conds, knots);
	//AccountConditions_1(A_1, b, conds, knots);
	vector<double> x(n_knots);
	LOS_LU(A);
	//LOC(A_1, b, x);

	knot k;
	k.x = 1;
	k.y = 0.25;

	rndPoint(k, elems[0], A.x);

	for (int i = 0; i < n_knots; i++)
		cout << FirstCondition(2, knots[i]) << "\t" << x[i] << "\t" <<
		abs(x[i] - FirstCondition(2, knots[i])) << endl;
	
	return 0;
}