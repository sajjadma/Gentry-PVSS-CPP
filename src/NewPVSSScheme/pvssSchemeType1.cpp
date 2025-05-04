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
    NTL::vec_ZZX operator*(const NTL::Mat<NTL::ZZX> &A, const NTL::vec_ZZX &b) {
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

    void _generateTrapdoor(NTL::mat_ZZ_pE &A, NTL::Mat<NTL::ZZX> &trapdoor, const long n, const long m,
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

        G = G - A_hat * NTL::conv<NTL::mat_ZZ_pE, NTL::Mat<NTL::ZZ_pX> >(
                NTL::conv<NTL::Mat<NTL::ZZ_pX>, NTL::Mat<NTL::ZZX> >(trapdoor));

        for (long i = 0; i < n; i++) {
            for (long j = 0; j < _n; j++) {
                A[i][_m + j] = G[i][j];
            }
        }
    }

    void _preSample(NTL::vec_ZZX &x, const NTL::Mat<NTL::ZZX> &trapdoor, const NTL::mat_ZZ_pE &A,
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

        params.k = securityParameter;
        power2(params.bound, params.k / 2 - 3);

        NTL::ZZ_p::init(params.q);

        clear(params.f);
        SetCoeff(params.f, params.k);
        SetCoeff(params.f, 0);
        NTL::ZZ_pE::init(params.f);

        clear(params.a);
        for (long i = 0; i < params.k; i++) {
            SetCoeff(params.a, i, RandomBnd(params.p));
        }

        params.w = params.threshold + 1 + 4 * params.numberOfParties;
        params.o = 3 * params.numberOfParties;
        params.v.SetLength(params.w);
        for (long i = 0; i < params.w; i++) {
            while (true) {
                random(params.v[i]);
                if (NTL::ZZ_pX unused; InvModStatus(unused, rep(params.v[i]), params.f) == 0) break;
            }
        }

        constexpr long n = 2;
        const long m = 2 * n * (securityParameter * 4 + 4);
        NTL::Mat<NTL::ZZX> td;
        _generateTrapdoor(params.A, td, n, m, params.q, params.f, params.bound);

        NTL::vec_ZZ_pE tmp;
        random(tmp, m);
        params.t = params.A * tmp;

        params.u.SetDims(params.w, params.w);
        for (long i = 0; i < params.w; i++) {
            for (long j = 0; j < params.w; j++) {
                if (i == j) continue;

                _preSample(params.u[i][j], td, params.A, params.v[i] / params.v[j] * params.t, params.q,
                           params.f, params.bound);
            }
        }

        NTL::ZZ_p::init(params.p);
        random(params.h, params.o);

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
        DistributionProof proof;
        NTL::vec_ZZX input;
        input.SetLength(params.w);
        NTL::Mat<NTL::ZZX> M;
        M.SetDims(params.o, params.w);
        long inputIndex = 0;
        long outputIndex = 0;
        proof.encryptedShares.SetLength(params.numberOfParties);

        // Polynomial part starts
        input[inputIndex] = secret;
        for (inputIndex = 1; inputIndex <= params.threshold; inputIndex++) {
            input[inputIndex] = to_ZZX(RandomBnd(params.p));
        }
        // Polynomial part ends

        for (long i = 0; i < params.numberOfParties; i++) {
            NTL::ZZ b = NTL::to_ZZ(1);
            NTL::ZZ share = NTL::to_ZZ(0);
            for (long j = 0; j < params.threshold + 1; j++) {
                M.put(outputIndex, j, to_ZZX(b));
                share += coeff(input[j], 0) * b;
                MulMod(b, b, i + 1, params.p);
            }

            input[inputIndex] = share;
            M.put(outputIndex, inputIndex, NTL::to_ZZX(-1));
            inputIndex++;
            outputIndex++;

            NTL::ZZX r, e1, e2;
            for (long j = 0; j < params.k + 1; j++) {
                SetCoeff(r, i, RandomBnd(params.bound));
                SetCoeff(e1, i, RandomBnd(params.bound));
                SetCoeff(e2, i, RandomBnd(params.bound));
            }
            input[inputIndex] = r;

            proof.encryptedShares.push_back(f1 * encodedShare + f2 * r);

            for (long j = 0; j < params.encryptionParams->cipherSize; j++) {
                outputIndex++;
                for (long k = 0; k < params.encryptionParams->plainSize; k++) {
                    const long index = params.threshold + 1 + (i * params.encryptionParams->plainSize) + k;
                    M1.put(outputIndex, index, f1[j][k]);
                }

                for (long k = 0; k < params.encryptionParams->randomSize; k++) {
                    for (long l = 0; l < params.vcParams->secondInputPartitionSize; l++) {
                        const long index = (i * params.encryptionParams->randomSize) +
                                           (k * params.vcParams->secondInputPartitionSize) + l;
                        M2.put(outputIndex, index,
                               l == 0
                                   ? f2[j][k]
                                   : f2[j][k] * params.vcParams->secondInputBound *
                                     M2[outputIndex][index - 1]);
                    }
                }
            }
        }

        VC::Auxiliary *auxiliary = new MyVectorCommitment::VectorCommitmentType2::MyAuxiliary();
        this->vectorCommitmentSystem->commit(proof.commitment, auxiliary, params.vcParams, firstInput, secondInput);
        this->vectorCommitmentSystem->open(proof.proof, params.vcParams, auxiliary, M1, M2);
    }

    bool verifyDistribution(const Params &params, const std::vector<PublicKey> &publicKeys,
                            const DistributionProof &proof) {
    }

    DecryptionProof decryptShare(const Params &params, const PublicKey &publicKey, const PrivateKey &privateKey,
                                 const Cipher &encryptedShare) {
    }

    bool verifyDecryption(const Params &params, const PublicKey &publicKey, const Cipher &encryptedShare,
                          const DecryptionProof &proof) {
        return true;
    }

    NTL::ZZ reconstruct(const Params &params, const std::vector<NTL::ZZ> &decryptedShares) {
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
