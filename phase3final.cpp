#include <bits/stdc++.h>
using namespace std;
const long double EPSILON = 1e-9;
int call = 0;
void printx(vector<long double> &x0)
{
    cout << "x0 = {";
    for (int i = 0; i < x0.size(); i++)
    {
        cout << x0[i] << (i == x0.size() - 1 ? "" : ", ");
    }
    cout << "}" << endl;
}

long double getRandomNumber(long double a, long double b, mt19937 &gen)
{
    // Check if a is greater than b, and swap if necessary
    if (a > b)
        swap(a, b);

    uniform_real_distribution<> dis(a, b); // Uniform distribution between a and b

    // Generate and return a random number between a and b
    return dis(gen);
}

long double mymin(long double a, long double b)
{
    if (a < b)
        return a;
    else
        return b;
}

long double mymax(long double a, long double b)
{
    if (a > b)
        return a;
    else
        return b;
}

int rankOfMatrix(vector<vector<long double>> &matrix)
{
    int rows = matrix.size();
    int cols = matrix[0].size();
    int rank = 0;

    for (int col = 0; col < cols; ++col)
    {
        // Find the pivot row
        int pivotRow = rank;
        while (pivotRow < rows && fabs(matrix[pivotRow][col]) < EPSILON)
        {
            ++pivotRow;
        }

        // If no pivot is found in this column, skip it
        if (pivotRow == rows)
            continue;

        // Swap the current row with the pivot row
        swap(matrix[rank], matrix[pivotRow]);

        // Normalize the pivot row
        long double pivotValue = matrix[rank][col];
        for (int j = col; j < cols; ++j)
        {
            matrix[rank][j] /= pivotValue;
        }

        // Eliminate all other rows using the pivot row
        for (int i = 0; i < rows; ++i)
        {
            if (i != rank && fabs(matrix[i][col]) > EPSILON)
            {
                long double factor = matrix[i][col];
                for (int j = col; j < cols; ++j)
                {
                    matrix[i][j] -= factor * matrix[rank][j];
                }
            }
        }

        // Move to the next row
        ++rank;
    }

    return rank;
}

// Function to check if a set of vectors is linearly dependent
bool isLinearlyDependent(vector<vector<long double>> &vectors)
{
    int numVectors = vectors[0].size();
    int rank = rankOfMatrix(vectors);
    // If the rank is less than the number of vectors, they are linearly dependent
    return rank < numVectors;
}

long double objective_func(int index, vector<long double> x)
{

    long double ans = 0;
    if (index == 1)
    {
        ans = (x[0] - 10) * (x[0] - 10) * (x[0] - 10) + (x[1] - 20) * (x[1] - 20) * (x[1] - 20);
    }
    else if (index == 2)
    {
        long double numerator = pow(sin(2 * M_PI * x[0]), 3) * sin(2 * M_PI * x[1]);
        long double denominator = pow(x[0], 3) * (x[0] + x[1]);
        ans = numerator / denominator;
    }
    else if (index == 3)
    {
        ans = x[0] + x[1] + x[2];
    }

    return ans;
}

long double penalty_term(int index, vector<long double> x)
{
    long double ans = 0;
    if (index == 1)
    {
        long double g1 = (x[0] - 5) * (x[0] - 5) + (x[1] - 5) * (x[1] - 5) - 100;
        long double g2 = (x[0] - 6) * (x[0] - 6) + (x[1] - 5) * (x[1] - 5) - 82.81;
        ans = mymin(0, g1) * mymin(0, g1) + mymax(0, g2) * mymax(0, g2);
    }
    else if (index == 2)
    {
        long double g1 = (x[0]) * (x[0]) - x[1] + 1;
        long double g2 = 1 - x[0] + (x[1] - 4) * (x[1] - 4);
        ans = mymax(0, g1) * mymax(0, g1) + mymax(0, g2) * mymax(0, g2);
    }
    else if (index == 3)
    {
        long double g1 = -1 + 0.0025 * (x[3] + x[5]);
        // g1 /= 4;

        long double g2 = -1 + 0.0025 * (-x[3] + x[4] + x[6]);
        // g2 /= 3.95;

        long double g3 = -1 + 0.01 * (-x[5] + x[7]);
        // g3 /= 8.9;

        long double g4 = 100 * x[0] - x[0] * x[5] + 833.33252 * x[3] - 83333.333;
        // g4 /= 1649999.187;

        long double g5 = x[1] * x[3] - x[1] * x[6] - 1250 * x[3] + 1250 * x[4];
        // g5 /= 9900000;

        long double g6 = x[2] * x[4] - x[2] * x[7] - 2500 * x[4] + 1250000;
        // g6 /= 8650000;

        ans = mymax(0, g1) * mymax(0, g1) + mymax(0, g2) * mymax(0, g2) + mymax(0, g3) * mymax(0, g3) +
              mymax(0, g4) * mymax(0, g4) + mymax(0, g5) * mymax(0, g5) + mymax(0, g6) * mymax(0, g6);
    }

    return ans;
}

long double penalty_func(int index, vector<long double> x, long double r)
{
    call++;
    return objective_func(index, x) + r * penalty_term(index, x);
}

vector<pair<long double, long double>> bounds(int index)
{

    vector<pair<long double, long double>> ans;
    if (index == 1)
    {
        pair<long double, long double> p1 = {13, 20};
        pair<long double, long double> p2 = {0, 4};
        ans.push_back(p1);
        ans.push_back(p2);
    }
    else if (index == 2)
    {
        pair<long double, long double> p1 = {0, 10};
        ans.push_back(p1);
        ans.push_back(p1);
    }
    else if (index == 3)
    {
        pair<long double, long double> p1 = {100, 10000};
        pair<long double, long double> p2 = {1000, 10000};
        pair<long double, long double> p3 = {10, 1000};
        ans.push_back(p1);
        ans.push_back(p2);
        ans.push_back(p2);
        ans.push_back(p3);
        ans.push_back(p3);
        ans.push_back(p3);
        ans.push_back(p3);
        ans.push_back(p3);
    }
    return ans;
}

long double norm(vector<long double> x)
{
    long double sum = 0;
    for (int i = 0; i < x.size(); i++)
    {
        sum += x[i] * x[i];
    }
    sum = sqrt(sum);
    return sum;
}

void normalize(vector<long double> &x)
{
    long double mag = norm(x);
    for (int i = 0; i < x.size(); i++)
    {
        x[i] /= mag;
    }
}

pair<vector<long double>, vector<long double>> boundingPhase(int index, vector<long double> x0, vector<long double> dir, long double r, long double delta = 0.15)
{

    int dim = x0.size();
    vector<long double> xleft(dim), xright(dim);
    for (int i = 0; i < x0.size(); i++)
    {
        xleft[i] = x0[i] - delta * dir[i];
        xright[i] = x0[i] + delta * dir[i];
    }

    long double f0 = penalty_func(index, x0, r);
    long double fleft = penalty_func(index, xleft, r);
    long double fright = penalty_func(index, xright, r);
    long double f1, f2;

    vector<long double> x1(dim), x2(dim);

    if (fleft <= f0 && f0 <= fright)
    {
        x1 = xleft;
        delta = -delta;
        f1 = fleft;
    }
    else if (fright <= f0 && f0 <= fleft)
    {
        x1 = xright;
        f1 = fright;
    }
    else if (f0 < fleft && f0 < fright)
    {
        return {xleft, xright};
    }
    // else
    // {
    //     x1 = xright;
    //     f1 = fright;
    // }

    for (int i = 0; i < dim; i++)
    {
        x2[i] = x1[i] + 2 * delta * dir[i];
    }

    f2 = penalty_func(index, x2, r);

    int k = 2;

    while (f0 >= f1 && f1 >= f2)
    {
        x0 = x1;
        x1 = x2;
        f0 = f1;
        f1 = f2;

        for (int i = 0; i < dim; i++)
        {
            x2[i] = x1[i] + pow(2, k) * delta * dir[i];
        }

        f2 = penalty_func(index, x2, r);
        k++;
    }

    return {x0, x2};
}

vector<long double> goldenSectionSearch(int index, vector<long double> a, vector<long double> b, long double r, long double eps = 1e-15)
{
    long double phi = 0.618;
    int dim = a.size();
    vector<long double> x1(dim), x2(dim);

    vector<long double> dir(dim);
    for (int i = 0; i < dim; i++)
    {
        dir[i] = b[i] - a[i];
    }

    for (int i = 0; i < dim; i++)
    {
        x1[i] = a[i] + phi * dir[i];
        x2[i] = b[i] - phi * dir[i];
    }

    long double f1 = penalty_func(index, x1, r);
    long double f2 = penalty_func(index, x2, r);

    while (norm(dir) >= eps)
    {
        if (f1 < f2)
        {
            a = x2;
            x2 = x1;
            f2 = f1;

            for (int i = 0; i < dim; i++)
            {
                dir[i] = b[i] - a[i];
            }

            for (int i = 0; i < dim; i++)
            {
                x1[i] = a[i] + phi * dir[i];
            }

            f1 = penalty_func(index, x1, r);
        }
        else
        {
            b = x1;
            x1 = x2;
            f1 = f2;

            for (int i = 0; i < dim; i++)
            {
                dir[i] = b[i] - a[i];
            }

            for (int i = 0; i < dim; i++)
            {
                x2[i] = b[i] - phi * dir[i];
            }

            f2 = penalty_func(index, x2, r);
        }
    }

    vector<long double> ans(dim);
    for (int i = 0; i < dim; i++)
    {
        ans[i] = (a[i] + b[i]) / 2;
    }
    return ans;
}

vector<long double> linsearch(int index, vector<long double> x0, vector<long double> dir, long double r)
{

    pair<vector<long double>, vector<long double>> p = boundingPhase(index, x0, dir, r);
    vector<long double> ans = goldenSectionSearch(index, p.first, p.second, r);
    return ans;
}

vector<long double> powellsDirectionMethod(int index, vector<long double> x0, long double r, long double eps = 1e-5, int maxIterations = 10000)
{

    int dim = x0.size();
    vector<vector<long double>> directions(dim, vector<long double>(dim, 0.0));
    for (int i = 0; i < dim - 1; i++)
    {
        directions[i][i + 1] = 1;
    }

    vector<long double> dir(dim, 0);
    dir[0] = 1;

    int k = 0;

    vector<long double> xnew = x0;

    while (norm(dir) > eps && k < maxIterations)
    {

        normalize(dir);

        for (int i = directions.size() - 1; i > 0; i--)
        {
            directions[i] = directions[i - 1];
        }
        directions[0] = dir;

        for (int i = 0; i < dim; i++)
        {

            xnew = linsearch(index, xnew, directions[i], r);
        }

        xnew = linsearch(index, xnew, directions[0], r);

        for (int i = 0; i < dim; i++)
        {
            dir[i] = xnew[i] - x0[i];
        }

        // Search for linear dependence
        x0 = xnew;
        k++;
        if (isLinearlyDependent(directions))
        {
            break;
        }
    }
    return xnew;
}

vector<long double> bracketOperatorMethod(int index, vector<long double> x0, long double r0, long double c = 12, long double eps = 1e-4, int maxIter = 1000)
{
    long double r = r0;
    vector<long double> xold = x0;
    vector<long double> xnew;
    long double pf, pi;
    int k = 0;
    // Open the file to write both iteration and function value
    ofstream file("convergence.csv");
    do
    {
        xnew = powellsDirectionMethod(index, xold, r);
        pi = penalty_func(index, xnew, r);
        r = r * c;
        xold = xnew;
        xnew = powellsDirectionMethod(index, xold, r);
        pf = penalty_func(index, xnew, r);
        xold = xnew;
        // printx(xnew);
        long double objectiveValue = objective_func(index, xnew);
        file << objectiveValue << "," << xnew[0] << "," << xnew[1] << endl; // Write iteration and function value
        k++;
    } while (abs(pf - pi) > eps && k < maxIter);
    cout << "Total number of iterations of Bracket Operator Method: " << k << endl;
    return xnew;
}

int main()
{
    int index, numvar;
    cout << "Enter the problem number: ";
    cin >> index;
    if (index == 1)
    {
        numvar = 2;
    }
    else if (index == 2)
    {
        numvar = 2;
    }
    else if (index == 3)
    {
        numvar = 8;
    }
    int runs = 10;
    int x = 10;
    // unsigned seed = std::chrono::steady_clock::now().time_since_epoch().count();
    unsigned seed = 42;
    mt19937 gen(seed);
    while (x--)
    {
        call=0;
        cout << "Run " << runs - x << ":" << endl;
        vector<long double> x0(numvar);
        vector<pair<long double, long double>> bound = bounds(index);
        for (int i = 0; i < numvar; i++)
        {
            x0[i] = getRandomNumber(bound[i].first, bound[i].second, gen);
        }
        cout << "Initial point is:" << endl;
        cout << "x0 = {";
        for (int i = 0; i < x0.size(); i++)
        {
            cout << x0[i] << (i == x0.size() - 1 ? "" : ", ");
        }
        cout << "}" << endl;
        vector<long double> xopt = bracketOperatorMethod(index, x0, 0.1);
        cout << "Optimized point is:" << endl;
        cout << "x_opt = {";
        for (int i = 0; i < xopt.size(); i++)
        {
            cout << xopt[i] << (i == xopt.size() - 1 ? "" : ", ");
        }
        cout << "}" << endl;

        // Print the function value
        cout << "f(x) = " << objective_func(index, xopt) << endl;
        cout << "Total number of Function Evaluations: " << call << endl;

        cout << endl;
        cout << endl;
    }

    return 0;
}
