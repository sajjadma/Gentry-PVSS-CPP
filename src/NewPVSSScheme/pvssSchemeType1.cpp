/* pvssSchemeType1.cpp
*
 * Copyright (C) 2025
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 **/

#include "pvssSchemeType1.hpp"

using namespace std;

namespace NewPVSSScheme::PVSSType1 {
    NTL::vec_ZZ g_inverse(const NTL::ZZ &u, const NTL::ZZ &q) {
        const long k = NumBits(q);
        NTL::vec_ZZ x, y;
        x.SetLength(k);
        y.SetLength(k);
        clear(x);
        clear(y);

        NTL::ZZ tmp = RandomBnd(q);
        x[k - 1] = tmp >= u ? NTL::ZZ(0) : NTL::ZZ(-1);

        NTL::ZZ u_curr = u;
        NTL::ZZ q_curr = q;

        for (long i = k - 2; i >= 0; i--) {
            u_curr -= bit(u, i + 1) << (i + 1);
            q_curr -= bit(q, i + 1) << (i + 1);

            NTL::ZZ c = -(u_curr + x[k - 1] * q);
            NTL::ZZ p, z;

            if (c < 0) {
                p = c + NTL::power2_ZZ(i + 1);
                z = -1;
            } else {
                p = c;
                z = 0;
            }

            tmp = RandomBnd(NTL::power2_ZZ(i + 1));
            x[i] = tmp < p ? (z + 1) : z;
        }

        for (long i = 1; i < k - 1; i++) {
            y[i] = 2 * x[i] - x[i - 1] + x[k - 1] * bit(q, i) + bit(u, i);
        }
        y[k - 1] = -x[k - 2] + x[k - 1] * bit(q, k - 1) + bit(u, k - 1);
        y[0] = 2 * x[0] + x[k - 1] * bit(q, 0) + bit(u, 0);

        return y;
    }

    void _generateTrapdoor(NTL::mat_ZZ_pE &A, NTL::mat_ZZ_pE &trapdoor, const long n, const long m, const NTL::ZZ &q,
                           const NTL::ZZ_pX &f, const NTL::ZZ &bound) {
        NTL::ZZ_pPush push;
        NTL::ZZ_p::init(q);
        NTL::ZZ_pE::init(f);

        const long k = NumBits(q);
        const long _n = n * k;
        const long _m = m - _n;

        NTL::mat_ZZ_pE G, A_hat;
        A_hat.SetDims(n, _m);
        G.SetDims(n, _n);
        A.SetDims(n, m);

        const NTL::ZZ a = SqrRoot(bound);
        NTL::ZZ_p::init(a);
        random(trapdoor, _m, _n);

        NTL::ZZ_p::init(q);
        for (long i = 0; i < n; i++) {
            for (long j = 0; j < _m; j++) {
                const NTL::ZZ_pE value = NTL::random_ZZ_pE();
                A_hat[i][j] = value;
                A[i][j] = value;
            }
        }

        for (long i = 0; i < n; i++) {
            for (long j = 0; j < k; j++) {
                const long index = i * k + j;
                G[i][index] = j == 0 ? NTL::to_ZZ_pE(1) : G[i][index - 1] * 2;
            }
        }

        G = G - A_hat * trapdoor;

        for (long i = 0; i < n; i++) {
            for (long j = 0; j < _n; j++) {
                A[i][_m + j] = G[i][j];
            }
        }
    }

    void _preSample(NTL::vec_ZZ_pE &x, const NTL::mat_ZZ_pE &trapdoor, const NTL::mat_ZZ_pE &A, const NTL::vec_ZZ_pE &b,
                    const NTL::ZZ &q, const NTL::ZZ_pX &f, const NTL::ZZ &bound) {
        NTL::ZZ_pPush push;
        NTL::ZZ_p::init(q);
        NTL::ZZ_pE::init(f);

        const long n = A.NumRows(), m = trapdoor.NumCols();
        const long k = m / n;
        const long d = deg(f);

        NTL::vec_ZZ_pE y;
        y.SetLength(m);
        const NTL::ZZ a = SqrRoot(bound);
        for (long i = 0; i < n; i++) {
            for (long j = 0; j < d; j++) {
                const NTL::ZZ u = rep(coeff(rep(b[i]), j));
                const NTL::vec_ZZ tmp = g_inverse(u, q);
                for (long l = 0; l < k; l++) {
                    const long index = i * k + l;
                    SetCoeff(y[index]._ZZ_pE__rep, j, to_ZZ_p(tmp[l]));
                }
            }
        }

        x = trapdoor * y;
        x.append(y);
    }

    Params setup(const long securityParameter, const long numberOfParties, const long threshold) {
        NTL::ZZ_pPush push;
        Params params;
        params.numberOfParties = numberOfParties;
        params.threshold = threshold;

        RandomPrime(params.p, securityParameter);
        RandomPrime(params.q, securityParameter * 4 + 4);

        params.k = securityParameter;
        power2(params.bound, params.k / 2 - 3);

        NTL::ZZ_p::init(params.q);

        SetCoeff(params.f, params.k);
        SetCoeff(params.f, 0);
        NTL::ZZ_pE::init(params.f);

        random(params.a);

        const long w = params.threshold + 1 + 4 * params.numberOfParties;
        const long o = 3 * params.numberOfParties;
        random(params.v, w);

        constexpr long n = 2;
        const long m = 2 * n * (securityParameter * 4 + 4);
        NTL::mat_ZZ_pE td;
        _generateTrapdoor(params.A, td, n, m, params.q, params.f, params.bound);

        NTL::vec_ZZ_pE tmp;
        random(tmp, w);
        params.t = params.A * tmp;

        params.u.SetDims(w, w);
        for (long i = 0; i < w; i++) {
            for (long j = 0; j < w; j++) {
                if (i == j) continue;

                _preSample(params.u[i][j], td, params.A, params.v[i] / params.v[j] * params.t, params.q,
                           params.f, params.bound);
            }
        }

        NTL::ZZ_p::init(params.p);
        random(params.h, o);

        return params;
    }

    KeyPair generateKey(const Params &params, long index) {
        NTL::ZZ_pPush push;
        KeyPair key;
        key.privateKey.a = params.a;
        key.publicKey.a = params.a;

        NTL::ZZ_p::init(params.bound);
        NTL::ZZ_pE::init(params.f);

        random(key.privateKey.s);

        NTL::ZZ_pE e;
        random(e);

        NTL::ZZ_p::init(params.q);

        key.publicKey.b = key.privateKey.a * key.privateKey.s + e;

        return key;
    }
}
