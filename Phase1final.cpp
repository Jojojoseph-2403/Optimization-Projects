// Roll No. : 210103012, 210103060
#include <bits/stdc++.h>
using namespace std;

int boundingPhaseIterations = 0; // Global counter for Bounding Phase Method iterations
int goldenSectionIterations = 0; // Global counter for Golden Section Search Method iterations

// Bounding Phase Method
pair<double, double> boundingPhaseMethod(function<double(double)> f, double x0, double a, double b, bool findMinimum = true, double delta = 0.001, int maxIter = 1000)
{
    srand(static_cast<unsigned int>(time(0))); // Seed for random number generation
    int k = 0;
    boundingPhaseIterations = 0; // Reset counter for each run

    double f_x0 = f(x0);
    double f_x0_minus = f(x0 - fabs(delta));
    double f_x0_plus = f(x0 + fabs(delta));

    // Determine the initial search direction
    if ((f_x0_minus >= f_x0 && f_x0 >= f_x0_plus) == findMinimum)
    {
        delta = fabs(delta);
    }
    else if ((f_x0_minus <= f_x0 && f_x0 <= f_x0_plus) == findMinimum)
    {
        delta = -fabs(delta);
    }
    else
    {
        return boundingPhaseMethod(f, x0, a, b, findMinimum, delta, maxIter); // Retry with a new random guess
    }

    double x_low = x0;
    double x_high = x0;

    // Iterative expansion
    while (k < maxIter)
    {
        ++boundingPhaseIterations; // Increment the counter
        double x_next = x0 + pow(2, k) * delta;

        if ((f(x_next) > f_x0) == findMinimum)
        {
            x_high = x_next;
            x_low = x0 - pow(2, k - 1) * delta;
            break;
        }

        x0 = x_next;
        f_x0 = f(x0);
        ++k;
    }


    if (k == maxIter)
    {
        cout << "Bounding Phase Method failed to find a suitable interval within max iterations." << endl;
        exit(1);
    }

    // Ensure the output range (x_low, x_high) satisfies x_low < x_high
    if (x_low > x_high)
    {
        swap(x_low, x_high);
    }

    // Return the bracketing interval
    return make_pair(x_low, x_high);
}

// Golden Section Search Method
pair<double, double> goldenSectionSearch(function<double(double)> f, double a, double b, bool findMinimum = true, double tolerance = 1e-5)
{
    goldenSectionIterations = 0;    // Reset counter for each run
    double phi = (sqrt(5) - 1) / 2; // Golden ratio

    double c = b - phi * (b - a);
    double d = a + phi * (b - a);

    while (abs(b - a) > tolerance)
    {
        ++goldenSectionIterations; // Increment the counter
        if ((findMinimum && f(c) < f(d)) || (!findMinimum && f(c) > f(d)))
        {
            b = d;
        }
        else
        {
            a = c;
        }

        c = b - phi * (b - a);
        d = a + phi * (b - a);
    }

    return make_pair(a, b); // Return the final interval
}

int main()
{
    int choice;
    cout << "Select a function to optimize (1-6):" << endl;
    cout << "1. f(x) = (2x-5)^4 - (x^2-1)^3" << endl;
    cout << "2. f(x) = 8 + x^3 - 2x - 2e^x" << endl;
    cout << "3. f(x) = 4x sin(x)" << endl;
    cout << "4. f(x) = 2(x-3)^2 + e^(0.5x^2)" << endl;
    cout << "5. f(x) = x^2 - 10e^(0.1x)" << endl;
    cout << "6. f(x) = 20 sin(x) - 15x^2" << endl;
    cin >> choice;

    function<double(double)> func;
    double a, b;
    bool findMinimum;
    string funcStr;

    switch (choice)
    {
    case 1:
        func = [](double x)
        { return pow(2 * x - 5, 4) - pow(x * x - 1, 3); };
        a = -10.0;
        b = 0.0;
        findMinimum = false;
        funcStr = "f(x) = (2x-5)^4 - (x^2-1)^3";
        break;
    case 2:
        func = [](double x)
        { return 8 + x * x * x - 2 * x - 2 * exp(x); };
        a = -2.0;
        b = 1.0;
        findMinimum = false;
        funcStr = "f(x) = 8 + x^3 - 2x - 2e^x";
        break;
    case 3:
        func = [](double x)
        { return 4 * x * sin(x); };
        a = 0.5;
        b = M_PI;
        findMinimum = false;
        funcStr = "f(x) = 4x sin(x)";
        break;
    case 4:
        func = [](double x)
        { return 2 * pow(x - 3, 2) + exp(0.5 * pow(x, 2)); };
        a = -2.0;
        b = 3.0;
        findMinimum = true;
        funcStr = "f(x) = 2(x-3)^2 + e^(0.5x^2)";
        break;
    case 5:
        func = [](double x)
        { return pow(x, 2) - 10 * exp(0.1 * x); };
        a = -6.0;
        b = 6.0;
        findMinimum = true;
        funcStr = "f(x) = x^2 - 10e^(0.1x)";
        break;
    case 6:
        func = [](double x)
        { return 20 * sin(x) - 15 * pow(x, 2); };
        a = -4.0;
        b = 4.0;
        findMinimum = false;
        funcStr = "f(x) = 20 sin(x) - 15x^2";
        break;
    default:
        cout << "Invalid choice." << endl;
        return 1;
    }

    double incrementX0 = (b - a) / 20; // Increment value for x0
    double delta = 0.001;              // Initial value for delta
    double x0 = a + 0.11111;

    // Run the selected function 10 times
    for (int i = 0; i < 10; ++i)
    {
        cout << "Run #" << (i + 1) << endl;

        // Output the function and its range
        cout << "Function: " << funcStr;
        cout << (findMinimum ? " (Minimize)" : " (Maximize)") << endl;
        cout << "Given range: (" << a << ", " << b << ")" << endl;

        // Output the initial values of x0 and delta
        cout << "Initial x0: " << x0 << endl;
        cout << "Initial delta: " << delta << endl;

        // Call bounding phase method
        auto boundingResult = boundingPhaseMethod(func, x0, a, b, findMinimum, delta);

        // Output the range obtained from bounding phase method
        cout << "Bounding Phase Method found the range: (" << boundingResult.first << ", " << boundingResult.second << ")" << endl;

        // Use the range from bounding phase method for Golden Section Search
        auto goldenResult = goldenSectionSearch(func, boundingResult.first, boundingResult.second, findMinimum);

        // Evaluate the function at the final point
        double optimalValue = func((goldenResult.first + goldenResult.second) / 2);

        // Output result from Golden Section Search
        cout << "Golden Section Search Method found the range for " << (findMinimum ? "minimum" : "maximum") << ": (" << goldenResult.first << ", " << goldenResult.second << ")" << endl;
        cout << "Optimal point: " << (goldenResult.first + goldenResult.second) / 2 << endl;
        cout << "Value of f(x) at optimal point: " << optimalValue << endl;

        // Output the number of iterations for each method
        cout << "Number of iterations in Bounding Phase Method: " << boundingPhaseIterations << endl;
        cout << "Number of iterations in Golden Section Search Method: " << goldenSectionIterations << endl;
        cout << "--------------------------------------" << endl;

        // Increment x0 and delta for the next run
        x0 += incrementX0;
        delta += 0.001;
    }

return 0;

}


