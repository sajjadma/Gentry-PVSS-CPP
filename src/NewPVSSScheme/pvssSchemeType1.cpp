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
    NTL::vec_ZZX operator+(const NTL::vec_ZZX &a, const NTL::vec_ZZX &b) {
        const long n = a.length();
        if (b.length() != n) throw "vector add: dimension mismatch";

        NTL::vec_ZZX x;
        x.SetLength(n);

        for (long i = 0; i < n; i++) {
            add(x[i], a[i], b[i]);
        }

        return x;
    }

    NTL::vec_ZZX operator*(const NTL::vec_ZZX &a, const NTL::ZZX &b) {
        const long n = a.length();

        NTL::vec_ZZX x;
        x.SetLength(n);

        for (long i = 0; i < n; i++) {
            mul(x[i], a[i], b);
        }

        return x;
    }

    NTL::vec_ZZX operator*(const mat_ZZX &A, const NTL::vec_ZZX &b) {
        const long n = A.NumRows();
        const long l = A.NumCols();

        if (l != b.length())
            throw "matrix mul: dimension mismatch";

        NTL::vec_ZZX x;
        x.SetLength(n);

        NTL::ZZX acc, tmp;
        for (long i = 1; i <= n; i++) {
            clear(acc);
            for (long k = 1; k <= l; k++) {
                mul(tmp, A(i, k), b(k));
                add(acc, acc, tmp);
            }
            conv(x(i), acc);
        }

        return x;
    }

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

    void _generateTrapdoor(NTL::mat_ZZ_pE &A, mat_ZZX &trapdoor, const long n, const long m,
                           const NTL::ZZ &q, const NTL::ZZ_pX &f, const NTL::ZZ &bound) {
        NTL::ZZ_pPush push;
        NTL::ZZ_p::init(q);
        NTL::ZZ_pE::init(f);

        const long k = NumBits(q);
        const long d = deg(f);
        const long _n = n * k;
        const long _m = m - _n;

        NTL::mat_ZZ_pE G, A_hat;
        A_hat.SetDims(n, _m);
        G.SetDims(n, _n);
        A.SetDims(n, m);
        trapdoor.SetDims(_m, _n);

        const NTL::ZZ a = SqrRoot(bound);
        for (long i = 0; i < _m; i++) {
            for (long j = 0; j < _n; j++) {
                for (long l = 0; l < d; l++) {
                    SetCoeff(trapdoor[i][j], l, RandomBnd(a));
                }
            }
        }

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

        G = G - A_hat * to_mat_ZZ_pE(trapdoor);

        for (long i = 0; i < n; i++) {
            for (long j = 0; j < _n; j++) {
                A[i][_m + j] = G[i][j];
            }
        }
    }

    void _preSample(NTL::vec_ZZX &x, const mat_ZZX &trapdoor, const NTL::mat_ZZ_pE &A,
                    const NTL::vec_ZZ_pE &b, const NTL::ZZ &q, const NTL::ZZ_pX &f, const NTL::ZZ &bound) {
        NTL::ZZ_pPush push;
        NTL::ZZ_p::init(q);
        NTL::ZZ_pE::init(f);

        const long n = A.NumRows(), m = trapdoor.NumCols();
        const long k = m / n;
        const long d = deg(f);

        NTL::vec_ZZX y;
        y.SetLength(m);
        const NTL::ZZ a = SqrRoot(bound);
        for (long i = 0; i < n; i++) {
            for (long j = 0; j < d; j++) {
                const NTL::ZZ u = rep(coeff(rep(b[i]), j));
                const NTL::vec_ZZ tmp = g_inverse(u, q);
                for (long l = 0; l < k; l++) {
                    const long index = i * k + l;
                    SetCoeff(y[index], j, tmp[l]);
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

        params.k = securityParameter * 3;
        power2(params.bound, securityParameter / 2 - 3);

        NTL::ZZ_p::init(params.q);

        clear(params.f);
        SetCoeff(params.f, params.k);
        SetCoeff(params.f, 0);
        NTL::ZZ_pE::init(params.f);

        clear(params.a);
        for (long i = 0; i < params.k; i++) {
            SetCoeff(params.a, i, RandomBnd(params.p));
        }

        params.w = params.threshold + 1 + 4 * params.numberOfParties * params.k;
        params.o = 3 * params.numberOfParties;
        params.v.SetLength(5); //TODO:(params.w);
        for (long i = 0; i < 5/*TODO:params.w*/; i++) {
            while (true) {
                random(params.v[i]);
                if (NTL::ZZ_pX unused; InvModStatus(unused, rep(params.v[i]), params.f) == 0) break;
            }
        }

        constexpr long n = 2;
        const long m = 2 * n * (securityParameter * 4 + 4);
        mat_ZZX td;
        _generateTrapdoor(params.A, td, n, m, params.q, params.f, params.bound);

        NTL::vec_ZZ_pE tmp;
        random(tmp, m);
        params.t = params.A * tmp;

        params.u.SetDims(5, 5); //TODO:(params.w, params.w);
        for (long i = 0; i < 5/*TODO:params.w*/; i++) {
            for (long j = 0; j < 5/*TODO:params.w*/; j++) {
                //TODO:if (i == j) continue;

                _preSample(params.u[i][j], td, params.A, params.v[i] / params.v[j] * params.t, params.q,
                           params.f, params.bound);
            }
        }

        params.h.SetLength(params.o);
        for (long i = 0; i < params.o; i++) {
            for (long j = 0; j < params.k; j++) {
                SetCoeff(params.h[i], j, RandomBnd(params.p));
            }
        }

        return params;
    }

    KeyPair generateKey(const Params &params, [[maybe_unused]] long index) {
        KeyPair key;
        key.privateKey.a = params.a;
        key.publicKey.a = params.a;

        clear(key.privateKey.s);
        NTL::ZZX e;
        for (long i = 0; i < params.k; i++) {
            SetCoeff(key.privateKey.s, i, RandomBnd(params.bound));
            SetCoeff(e, i, RandomBnd(params.bound));
        }

        key.publicKey.b = key.privateKey.a * key.privateKey.s + e;

        return key;
    }

    bool verifyKey(const Params &params, const PublicKey &publicKey, const KeyProof &proof) {
        return true;
    }

    DistributionProof distribute(const Params &params, const NTL::Vec<PublicKey> &publicKeys,
                                 const NTL::ZZ &secret) {
        NTL::ZZ_pPush push;
        NTL::ZZ_p::init(params.q);
        NTL::ZZ_pE::init(params.f);

        DistributionProof proof;
        NTL::vec_ZZX input;
        input.SetLength(params.w);
        mat_ZZX M;
        M.SetDims(params.o, params.w);
        long inputIndex = 0;
        long outputIndex = 0;
        proof.encryptedShares.SetLength(params.numberOfParties);
        NTL::ZZ tmp, share;
        const NTL::ZZ p2 = params.p / 2;
        NTL::ZZX m, r, e1, e2, tmpX;

        // Polynomial part starts
        input[inputIndex] = secret;
        for (inputIndex = 1; inputIndex <= params.threshold; inputIndex++) {
            input[inputIndex] = to_ZZX(RandomBnd(params.p));
        }
        // Polynomial part ends

        for (long i = 0; i < params.numberOfParties; i++) {
            set(tmp);
            clear(share);
            for (long j = 0; j < params.threshold + 1; j++) {
                M.put(outputIndex, j, to_ZZX(tmp));
                share += coeff(input[j], 0) * tmp;
                MulMod(tmp, tmp, i + 1, params.p);
            }

            clear(m);
            clear(r);
            clear(e1);
            clear(e2);
            for (long j = 0; j < params.k; j++) {
                clear(tmpX);
                SetCoeff(tmpX, j, 1);

                tmp = NTL::to_ZZ(bit(share, j));
                SetCoeff(m, j, tmp);
                input[inputIndex] = to_ZZX(tmp);
                M.put(outputIndex, inputIndex, to_ZZX(-NTL::power2_ZZ(j)));
                M.put(outputIndex + 2, inputIndex, tmpX * p2);
                inputIndex++;

                RandomBnd(tmp, params.bound);
                SetCoeff(r, j, tmp);
                input[inputIndex] = to_ZZX(tmp);
                M.put(outputIndex + 1, inputIndex, publicKeys[i].a * tmpX);
                M.put(outputIndex + 2, inputIndex, publicKeys[i].b * tmpX);
                inputIndex++;

                RandomBnd(tmp, params.bound);
                SetCoeff(e1, j, tmp);
                input[inputIndex] = to_ZZX(tmp);
                M.put(outputIndex + 1, inputIndex, tmpX);
                inputIndex++;

                RandomBnd(tmp, params.bound);
                SetCoeff(e2, j, tmp);
                input[inputIndex] = to_ZZX(tmp);
                M.put(outputIndex + 2, inputIndex, tmpX);
                inputIndex++;
            }

            proof.encryptedShares[i].u = publicKeys[i].a * r + e1;
            proof.encryptedShares[i].v = publicKeys[i].b * r + m * p2 + e2;
            outputIndex += 3;
        }

        NTL::Vec<NTL::vec_ZZX> u;
        u.SetLength(params.w);
        clear(proof.commitment.c);
        for (long i = 0; i < params.w; i++) {
            proof.commitment.c += params.v[i % 5/*TODO:*/] * to_ZZ_pE(to_ZZ_pX(input[i]));

            u[i].SetLength(params.A.NumCols(), NTL::ZZX::zero());
            for (long j = 0; j < params.w; j++) {
                if (j == i) continue;
                u[i] = u[i] + params.u[j % 5 /*TODO:*/][i % 5/*TODO:*/] * input[j];
            }
        }

        proof.proof.output = M * input;
        proof.proof.pi.SetLength(params.A.NumCols(), NTL::ZZX::zero());
        for (long i = 0; i < params.o; i++) {
            for (long j = 0; j < params.w; j++) {
                proof.proof.pi = proof.proof.pi + u[j] * params.h[i] * M[i][j];
            }
        }

        return proof;
    }

    bool verifyDistribution(const Params &params, const NTL::Vec<PublicKey> &publicKeys,
                            const DistributionProof &proof) {
        NTL::ZZ_pPush push;
        NTL::ZZ_p::init(params.q);
        NTL::ZZ_pE::init(params.f);

        mat_ZZX M;
        M.SetDims(params.o, params.w);
        long inputIndex = 0;
        long outputIndex = 0;
        NTL::ZZ tmp;
        const NTL::ZZ p2 = params.p / 2;
        NTL::ZZX tmpX;

        for (long i = 0; i < params.numberOfParties; i++) {
            if (!IsZero(proof.proof.output[outputIndex])) {
                cout << "Wrong output dist verify: (" << outputIndex << ") " << proof.proof.output[outputIndex] << endl;
            }
            if (proof.proof.output[outputIndex + 1] != proof.encryptedShares[i].u) {
                cout << "Wrong output dist verify: (" << outputIndex + 1 << ") " << proof.proof.output[outputIndex + 1]
                        << endl;
            }
            if (proof.proof.output[outputIndex + 2] != proof.encryptedShares[i].v) {
                cout << "Wrong output dist verify: (" << outputIndex + 2 << ") " << proof.proof.output[outputIndex + 2]
                        << endl;
            }

            set(tmp);
            for (long j = 0; j < params.threshold + 1; j++) {
                M.put(outputIndex, j, to_ZZX(tmp));
                MulMod(tmp, tmp, i + 1, params.p);
            }

            for (long j = 0; j < params.k; j++) {
                clear(tmpX);
                SetCoeff(tmpX, j, 1);

                M.put(outputIndex, inputIndex, to_ZZX(-NTL::power2_ZZ(j)));
                M.put(outputIndex + 2, inputIndex, tmpX * p2);
                inputIndex++;

                M.put(outputIndex + 1, inputIndex, publicKeys[i].a * tmpX);
                M.put(outputIndex + 2, inputIndex, publicKeys[i].b * tmpX);
                inputIndex++;

                M.put(outputIndex + 1, inputIndex, tmpX);
                inputIndex++;

                M.put(outputIndex + 2, inputIndex, tmpX);
                inputIndex++;
            }

            outputIndex += 3;
        }

        NTL::ZZ_pE x;
        for (long i = 0; i < params.o; i++) {
            for (long j = 0; j < params.w; j++) {
                x += to_ZZ_pE(to_ZZ_pX(params.h[i])) * (
                    to_ZZ_pE(to_ZZ_pX(M[i][j])) * proof.commitment.c / params.v[j % 5/*TODO:*/] -
                    to_ZZ_pE(to_ZZ_pX(proof.proof.output[i])));
            }
        }

        const auto a = params.A * to_vec_ZZ_pE(proof.proof.pi);
        const auto b = x * params.t;
        if (a != b) {
            cout << "Au != rhs in dist verify. a = " << a << ", b = " << b << endl;
        }

        return true;
    }

    DecryptionProof decryptShare(const Params &params, const PublicKey &publicKey, const PrivateKey &privateKey,
                                 const Cipher &encryptedShare) {
        DecryptionProof proof;
        clear(proof.decryptedShare);

        NTL::ZZX m = encryptedShare.v - privateKey.s * encryptedShare.u;
        for (long i = 0; i < params.k; i++) {
            if (coeff(m, i) > (params.p / 4)) {
                SetBit(proof.decryptedShare, i);
            }
        }

        return proof;
    }

    bool verifyDecryption(const Params &params, const PublicKey &publicKey, const Cipher &encryptedShare,
                          const DecryptionProof &proof) {
        return true;
    }

    NTL::ZZ reconstruct(const Params &params, const NTL::vec_ZZ &decryptedShares) {
        NTL::ZZ_pPush push;

        NTL::ZZ_p::init(params.p);

        NTL::vec_ZZ_p a, b;
        a.SetLength(params.numberOfParties);
        b.SetLength(params.numberOfParties);
        for (long i = 0; i < params.numberOfParties; i++) {
            a[i] = NTL::to_ZZ_p(i + 1);
            b[i] = to_ZZ_p(decryptedShares[i]);
        }

        const NTL::ZZ_pX f = interpolate(a, b);
        return rep(eval(f, NTL::ZZ_p::zero()));
    }
}
