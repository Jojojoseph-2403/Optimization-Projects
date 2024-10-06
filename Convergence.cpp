#include <bits/stdc++.h>
using namespace std;

const double EPSILON = 1e-9; // Tolerance for comparing floating-point numbers

// Function to perform Gaussian elimination and calculate the rank of the matrix
int rankOfMatrix(vector<vector<double>> &matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size();
    int rank = 0;

    for (int col = 0; col < cols; ++col) {
        int pivotRow = rank;
        while (pivotRow < rows && fabs(matrix[pivotRow][col]) < EPSILON) {
            ++pivotRow;
        }

        if (pivotRow == rows) continue;

        swap(matrix[rank], matrix[pivotRow]);
        double pivotValue = matrix[rank][col];

        for (int j = col; j < cols; ++j) {
            matrix[rank][j] /= pivotValue;
        }

        for (int i = 0; i < rows; ++i) {
            if (i != rank && fabs(matrix[i][col]) > EPSILON) {
                double factor = matrix[i][col];
                for (int j = col; j < cols; ++j) {
                    matrix[i][j] -= factor * matrix[rank][j];
                }
            }
        }

        ++rank;
    }

    return rank;
}

// Function to check if a set of vectors is linearly dependent
bool isLinearlyDependent(vector<vector<double>> &vectors) {
    int numVectors = vectors[0].size();
    int rank = rankOfMatrix(vectors);
    return rank < numVectors;
}

// Objective function calculations
double objective_func(int index, vector<double> x) {
    double ans = 0;
    int dim = x.size();
    if (index == 1) {
        for (int i = 0; i < dim; i++) {
            ans += (i + 1) * (x[i] * x[i]);
        }
    } else if (index == 2) {
        for (int i = 0; i < dim - 1; i++) {
            ans += 100 * pow(x[i + 1] - pow(x[i], 2), 2) + pow(x[i] - 1, 2);
        }
    } else if (index == 3) {
        for (int i = 1; i < dim; i++) {
            ans += (i + 1) * pow(2 * pow(x[i], 2) - x[i - 1], 2);
        }
        ans += pow(x[0] - 1, 2);
    } else if (index == 4) {
        for (int i = 0; i < dim; i++) {
            ans += pow(x[i] - 1, 2);
            if (i > 0) {
                ans -= x[i] * x[i - 1];
            }
        }
    } else if (index == 5) {
        double sum = 0;
        for (int i = 0; i < dim; i++) {
            ans += pow(x[i], 2);
            sum += 0.5 * (i + 1) * x[i];
        }
        ans += pow(sum, 2) + pow(sum, 4);
    }
    return ans;
}

// Bounds for each function
pair<double, double> bounds(int index, int n) {
    double ulimit, llimit;
    if (index == 1) {
        llimit = -5.12;
        ulimit = 5.12;
    } else if (index == 2) {
        llimit = -2.048;
        ulimit = 2.048;
    } else if (index == 3) {
        llimit = -10;
        ulimit = 10;
    } else if (index == 4) {
        llimit = -n * n;
        ulimit = -llimit;
    } else {
        llimit = -5;
        ulimit = 10;
    }
    return {llimit, ulimit};
}

// Calculate the norm of a vector
double norm(vector<double> x) {
    double sum = 0;
    for (int i = 0; i < x.size(); i++) {
        sum += x[i] * x[i];
    }
    return sqrt(sum);
}

// Normalize a vector
void normalize(vector<double> &x) {
    double mag = norm(x);
    for (int i = 0; i < x.size(); i++) {
        x[i] /= mag;
    }
}

// Bounding phase function
pair<vector<double>, vector<double>> boundingPhase(int index, vector<double> x0, vector<double> dir, double delta = 0.15) {
    int dim = x0.size();
    vector<double> xleft(dim), xright(dim);
    for (int i = 0; i < dim; i++) {
        xleft[i] = x0[i] - delta * dir[i];
        xright[i] = x0[i] + delta * dir[i];
    }

    double f0 = objective_func(index, x0);
    double fleft = objective_func(index, xleft);
    double fright = objective_func(index, xright);
    double f1, f2;
    vector<double> x1, x2(dim);

    if (fleft <= f0 && f0 <= fright) {
        x1 = xleft;
        delta = -delta;
        f1 = fleft;
    } else if (fright <= f0 && f0 <= fleft) {
        x1 = xright;
        f1 = fright;
    } else if (f0 < fleft && f0 < fright) {
        return {xleft, xright};
    } else {
        x1 = xright;
        f1 = fright;
    }

    for (int i = 0; i < dim; i++) {
        x2[i] = x1[i] + 2 * delta * dir[i];
    }

    f2 = objective_func(index, x2);
    int k = 2;

    while (f0 >= f1 && f1 >= f2) {
        x0 = x1;
        x1 = x2;
        f0 = f1;
        f1 = f2;

        for (int i = 0; i < dim; i++) {
            x2[i] = x1[i] + pow(2, k) * delta * dir[i];
        }

        f2 = objective_func(index, x2);
    }

    return {x0, x2};
}

// Golden section search
vector<double> goldenSectionSearch(int index, vector<double> a, vector<double> b, double eps = 1e-15) {
    double phi = 0.618;
    int dim = a.size();
    vector<double> x1(dim), x2(dim);
    vector<double> dir(dim);
    
    for (int i = 0; i < dim; i++) {
        dir[i] = b[i] - a[i];
    }

    for (int i = 0; i < dim; i++) {
        x1[i] = a[i] + phi * dir[i];
        x2[i] = b[i] - phi * dir[i];
    }

    double f1 = objective_func(index, x1);
    double f2 = objective_func(index, x2);

    while (norm(dir) >= eps) {
        if (f1 < f2) {
            a = x2;
            x2 = x1;
            f2 = f1;

            for (int i = 0; i < dim; i++) {
                dir[i] = b[i] - a[i];
            }

            for (int i = 0; i < dim; i++) {
                x1[i] = a[i] + phi * dir[i];
            }

            f1 = objective_func(index, x1);
        } else {
            b = x1;
            x1 = x2;
            f1 = f2;

            for (int i = 0; i < dim; i++) {
                dir[i] = b[i] - a[i];
            }

            for (int i = 0; i < dim; i++) {
                x2[i] = b[i] - phi * dir[i];
            }

            f2 = objective_func(index, x2);
        }
    }

    vector<double> ans(dim);
    for (int i = 0; i < dim; i++) {
        ans[i] = (a[i] + b[i]) / 2;
    }
    return ans;
}

// Line search using bounding phase and golden section search
vector<double> linsearch(int index, vector<double> x0, vector<double> dir, double delta = 0.12, double eps = 1e-7) {
    pair<vector<double>, vector<double>> p = boundingPhase(index, x0, dir, delta);
    vector<double> ans = goldenSectionSearch(index, p.first, p.second, eps);
    return ans;
}

// Powell's Direction Method
vector<double> powellsDirectionMethod(int index, vector<double> x0, double eps = 1e-15, int maxIterations = 1000) {
    int dim = x0.size();
    vector<vector<double>> directions(dim, vector<double>(dim, 0.0));
    for (int i = 0; i < dim - 1; i++) {
        directions[i][i + 1] = 1;
    }

    vector<double> dir(dim, 0);
    dir[0] = 1;

    int k = 0;
    vector<double> xnew = x0;
    vector<double> history; // To store objective function values for plotting

    // Open the file to write both iteration and function value
    ofstream file("convergence.csv");
    
    while (norm(dir) > eps && k < maxIterations) {
        normalize(dir);

        for (int i = directions.size() - 1; i > 0; i--) {
            directions[i] = directions[i - 1];
        }
        directions[0] = dir;

        for (int i = 0; i < dim; i++) {
            xnew = linsearch(index, xnew, directions[i]);
        }

        xnew = linsearch(index, xnew, directions[0]);

        for (int i = 0; i < dim; i++) {
            dir[i] = xnew[i] - x0[i];
        }

        // Record objective function value and iteration number
        double objectiveValue = objective_func(index, xnew);
        file << k << "," << objectiveValue << endl; // Write iteration and function value

        x0 = xnew;
        k++;
        if (isLinearlyDependent(directions)) {
            break;
        }
    }

    file.close(); // Close the file
    return xnew;
}

int main() {
    int index;
    double lowerBound, upperBound;
    int numvariables;

    // Ask the user for the function flag
    cout << "Select the function:" << endl;
    cout << "1. Sum of squares" << endl;
    cout << "2. Rosenbrock" << endl;
    cout << "3. Dixon-Price" << endl;
    cout << "4. Trid" << endl;
    cout << "5. Zakharov" << endl;
    cin >> index;

    // Ask the user for the number of variables
    cout << "Enter number of variables: ";
    cin >> numvariables;

    pair<double, double> p = bounds(index, numvariables);
    lowerBound = p.first;
    upperBound = p.second;

    // Create the initial guess vector based on the min and max value
    double xinitial = lowerBound + (upperBound - lowerBound) / 3.3;
    vector<double> x0(numvariables, xinitial);
    
    cout << "Initial point is:" << endl;
    cout << "x0 = {";
    for (int i = 0; i < x0.size(); i++) {
        cout << x0[i] << (i == x0.size() - 1 ? "" : ", ");
    }
    cout << "}" << endl;

    // Call Powell's Direction Method
    vector<double> result = powellsDirectionMethod(index, x0);

    // Print the result
    cout << "Optimized point is:" << endl;
    cout << "x_opt = {";
    for (int i = 0; i < result.size(); i++) {
        cout << result[i] << (i == result.size() - 1 ? "" : ", ");
    }
    cout << "}" << endl;

    // Print the function value
    cout << "f(x) = " << objective_func(index, result) << endl;
    return 0;
}
