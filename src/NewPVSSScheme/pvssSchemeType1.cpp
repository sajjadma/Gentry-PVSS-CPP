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
    Params setup(const long securityParameter, const long numberOfParties, const long threshold) {
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

        long n = 2;
        long m = 2 * n * (securityParameter * 4 + 4);
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
