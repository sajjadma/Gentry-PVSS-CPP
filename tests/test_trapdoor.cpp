#include <iostream>
#include <NTL/ZZ.h>
#include <NTL/ZZ_pXFactoring.h>

#define DEBUGGING

#include "tests.hpp"
#include "utils.hpp"
#include "pvssSchemeType1.hpp"

using namespace std;
using namespace NewPVSSScheme::PVSSType1;
using namespace NTL;

bool test_trapdoor_generation() {
    long securityParameter = 64;
    long k = securityParameter;
    long n = 2;
    long m = 2 * n * (securityParameter * 4 + 4);
    ZZ q;
    ZZ_pX f;
    ZZ bound;

    // Initialize parameters
    RandomPrime(q, securityParameter * 4 + 4);
    power2(bound, k / 2 - 3);
    ZZ_p::init(q);
    SetCoeff(f, k);
    SetCoeff(f, 0);
    ZZ_pE::init(f);

    mat_ZZ_pE A;
    Mat<ZZX> trapdoor;

    _generateTrapdoor(A, trapdoor, n, m, q, f, bound);

    // Verify matrix dimensions
    if (A.NumRows() != n || A.NumCols() != m) {
        return false;
    }
    if (trapdoor.NumRows() != m - n*NumBits(q) || trapdoor.NumCols() != n*NumBits(q)) {
        return false;
    }

    // Verify trapdoor entries are bounded
    long d = deg(f);
    for (long i = 0; i < trapdoor.NumRows(); i++) {
        for (long j = 0; j < trapdoor.NumCols(); j++) {
            for (long k = 0; k <= deg(trapdoor[i][j]); k++) {
                if (coeff(trapdoor[i][j], k) >= SqrRoot(bound)) {
                    return false;
                }
            }
        }
    }

    return true;
}

bool test_presample() {
    long securityParameter = 64;
    long k = securityParameter;
    long n = 2;
    long m = 2 * n * (securityParameter * 4 + 4);
    ZZ q;
    ZZ_pX f;
    ZZ bound;

    // Initialize parameters
    RandomPrime(q, securityParameter * 4 + 4);
    power2(bound, k / 2 - 3);
    ZZ_p::init(q);
    SetCoeff(f, k);
    SetCoeff(f, 0);
    ZZ_pE::init(f);

    mat_ZZ_pE A;
    Mat<ZZX> trapdoor;
    _generateTrapdoor(A, trapdoor, n, m, q, f, bound);

    vec_ZZ_pE b;
    random(b, n);

    vec_ZZX x;
    _preSample(x, trapdoor, A, b, q, f, bound);

    // Verify output dimensions
    if (x.length() != m) {
        return false;
    }

    // Convert x to ZZ_pE and verify A*x ≈ b
    ZZ_pPush push;
    ZZ_p::init(q);
    ZZ_pE::init(f);

    vec_ZZ_pE Ax;
    Ax.SetLength(n);
    for (long i = 0; i < n; i++) {
        for (long j = 0; j < m; j++) {
            ZZ_pX xj_px;
            conv(xj_px, x[j]);
            ZZ_pE xj;
            conv(xj, xj_px);
            Ax[i] += A[i][j] * xj;
        }
    }

    // Check Ax = b (mod q)
    for (long i = 0; i < n; i++) {
        ZZ_pE diff = Ax[i] - b[i];
        if (deg(rep(diff)) != -1) {
            return false;
        }
    }

    return true;
}

int main(int, char**) {
    if (!test_trapdoor_generation() || !test_presample()) {
        cout << LWEVSS_TESTS::failed << endl;
        return 1;
    }
    cout << LWEVSS_TESTS::passed << endl;
    return 0;
}