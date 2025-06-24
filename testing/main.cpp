#include <fstream>
#include <iostream>
using namespace std;

int main() {
    ofstream file("C:/Users/frodo/Desktop/test_output.csv");
    if (!file.is_open()) {
        cerr << "Failed to open file.\n";
        return 1;
    }
    file << "hello,world\n";
    file.close();
    cout << "Wrote file.\n";
    return 0;
}