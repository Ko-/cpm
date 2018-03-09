#include <iostream>
#include <cstdint>
using namespace std;

typedef uint16_t UINT16;
typedef uint32_t UINT32;

int canonical(UINT16 vec)
{
    UINT16 rolling = vec;
    for (int i = 1; i < 16; i++) {
        rolling = (rolling << 1) ^ (rolling >> 15);
        if (rolling < vec) return 0;
    }
    return 1;

}

//help tool for choosing circular Z
void testCPM()
{
    UINT16 pattern = 0x13;
    UINT16 histogram[17*16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    UINT16 nrtotest = 0;
    cout << "{";
    for (UINT32 i = 1; i < 0x10000; i++) { //or all parities
        UINT16 HWin = 0;
        if (canonical(i)) { //only consider the version with lowest z-coordinate
            UINT16 outcome = 0;
            UINT16 workI = i;
            for (UINT16 c = 1; workI > 0; c++, workI >>= 1) {
                if ((workI & 1) > 0) {
                    outcome ^= (pattern << c) ^ (pattern >> (16-c));
                    HWin += 1; //hw of parity, nr odd columns
                }
            }
            UINT16 HWout = 0; //hw of outcome, nr affected columns
            UINT16 workOutcome = outcome;
            while (workOutcome > 0) {
                if ((workOutcome & 1) > 0) HWout += 1;
                workOutcome >>= 1;
            }
            histogram[(HWin-1) + 16*HWout] += 1;
            UINT16 weight = 4*HWout;
            workOutcome = i & ~outcome;
            while(workOutcome > 0) {
                if ((workOutcome & 1) > 0) weight += 2;
                workOutcome >>= 1;
            }
            if (weight <= 25) { //maxweight
                cout << "{" << i << "," << outcome << "},";
                nrtotest += 1;
            }
        }
    }
    cout << "}" << endl << "To test: " << nrtotest << endl;
    cout << "Pattern: 0x" << hex << pattern << dec << endl;
    for (UINT16 HWout=0 ; HWout<17 ; HWout++ ){
        cout << "HW out " << HWout << "  ";
            for (UINT16 HWin = 1; HWin < 17; HWin++) {
                cout << " " << histogram[(HWin - 1) + 16 * HWout];
            }
        cout << endl;
    }
}

int main(void) {
    testCPM();
}
