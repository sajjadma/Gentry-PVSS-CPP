/* my_vector_commitment_type2.cpp
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

#include <my_implementation.hpp>
#include <sodium/crypto_generichash.h>

using namespace std;

namespace MyVectorCommitment {
    void VectorCommitmentType2::_generateTrapdoor(NTL::mat_ZZ &A, NTL::mat_ZZ &trapdoor, const long n, const long m,
                                                  const NTL::ZZ &module, const NTL::ZZ &bound) {
        const long k = NTL::NumBits(module);
        const long _n = n * k;
        const long _m = m - _n;

        NTL::mat_ZZ G, A_hat;
        A_hat.SetDims(n, _m);
        G.SetDims(n, _n);
        A.SetDims(n, m);
        trapdoor.SetDims(_m, _n);

        const NTL::ZZ a = SqrRoot(bound);
        for (long i = 0; i < _m; i++) {
            for (long j = 0; j < _n; j++) {
                trapdoor[i][j] = NTL::RandomBnd(a);
            }
        }

        for (long i = 0; i < n; i++) {
            for (long j = 0; j < _m; j++) {
                const NTL::ZZ value = NTL::RandomBnd(module);
                A_hat[i][j] = value;
                A[i][j] = value;
            }
        }

        for (long i = 0; i < n; i++) {
            for (long j = 0; j < k; j++) {
                const long index = i * k + j;
                G[i][index] = j == 0 ? NTL::to_ZZ(1) : G[i][index - 1] * 2;
            }
        }

        G = G - A_hat * trapdoor;

        for (long i = 0; i < n; i++) {
            for (long j = 0; j < _n; j++) {
                A[i][_m + j] = G[i][j];
            }
        }
    }

    void VectorCommitmentType2::_preSample(NTL::vec_ZZ &x, const NTL::mat_ZZ &trapdoor, const NTL::mat_ZZ &A,
                                           const NTL::vec_ZZ &b, const NTL::ZZ &bound) {
        const long n = A.NumRows(), m = trapdoor.NumCols();
        const long k = m / n;

        NTL::vec_ZZ y;
        y.SetLength(m);
        const NTL::ZZ a = NTL::SqrRoot(bound);
        for (long i = 0; i < n; i++) {
            NTL::ZZ u = b[i];
            for (long j = 0; j < k; j++) {
                const long index = i * k + j;
                y[index] = NTL::RandomBnd(a);
                if (u % 2 != y[index] % 2) {
                    y[index] += 1;
                }

                u = (u - y[index]) / 2;
            }
        }

        x = trapdoor * y;
        x.append(y);
    }

    // TODO: Set Parameters
    void VectorCommitmentType2::setup(const MyFramework::VC::Params *params, const long securityParameter,
                                      const long firstInputSize, const long secondInputSize, const long outputSize,
                                      const NTL::ZZ firstInputBound, const NTL::ZZ secondInputBound,
                                      const NTL::ZZ coefficientBound) {
        auto start = chrono::steady_clock::now();
        cout << "VC Setup starts" << endl;

        const auto myParams = (MyParams *) params;
        const long n = NTL::NumBits(secondInputBound);
        const long k = NTL::NumBits(
            (firstInputBound > coefficientBound ? firstInputBound : coefficientBound) * firstInputSize * outputSize);
        myParams->firstInputSize = firstInputSize;
        myParams->secondInputSize = secondInputSize * n;
        myParams->secondInputPartitionSize = n;
        myParams->outputSize = outputSize;
        myParams->firstInputBound = NTL::power2_ZZ(k);
        myParams->secondInputBound = 2;
        myParams->coefficientBound = myParams->firstInputBound;
        NTL::vec_ZZ v1Inv, v2Inv;
        v1Inv.SetLength(myParams->firstInputSize);
        myParams->v1.SetLength(myParams->firstInputSize);
        v2Inv.SetLength(myParams->secondInputSize);
        myParams->v2.SetLength(myParams->secondInputSize);
        myParams->h.SetLength(myParams->outputSize);
        myParams->h2.SetLength(myParams->secondInputSize);
        myParams->u.SetDims(myParams->firstInputSize + 3,//myParams->secondInputSize,
                            myParams->firstInputSize + 3);//myParams->secondInputSize);
        myParams->n = 8;
        myParams->l = 32 * myParams->n * k;
        myParams->p = NTL::power2_ZZ(4 * k);
        myParams->q = NTL::power2_ZZ(16 * k);
        myParams->beta = NTL::power2_ZZ(k);
        myParams->alpha = NTL::power2_ZZ(12 * k);
        myParams->t.SetLength(myParams->n);

        for (long i = 0; i < myParams->firstInputSize; i++) {
            if (i > 2) {
                v1Inv[i] = v1Inv[i % 3];
                myParams->v1[i] = myParams->v1[i % 3];
                continue;
            }
            NTL::ZZ tmp;
            do {
                tmp = NTL::RandomBnd(myParams->q);
            } while (NTL::GCD(tmp, myParams->q) != 1);

            v1Inv[i] = tmp;
            NTL::InvMod(myParams->v1[i], tmp, myParams->q);
        }

        for (long i = 0; i < myParams->secondInputSize; i++) {
            if (i > 2) {
                v2Inv[i] = v2Inv[i % 3];
                myParams->v2[i] = myParams->v2[i % 3];
                continue;
            }
            NTL::ZZ tmp;
            do {
                tmp = NTL::RandomBnd(myParams->p);
            } while (NTL::GCD(tmp, myParams->q) != 1);

            v2Inv[i] = tmp;
            NTL::InvMod(myParams->v2[i], tmp, myParams->q);
        }

        for (long i = 0; i < myParams->outputSize; i++) {
            myParams->h[i] = NTL::RandomBnd(myParams->p);
        }

        for (long i = 0; i < myParams->secondInputSize; i++) {
            myParams->h2[i] = NTL::RandomBnd(myParams->p);
        }

        for (long i = 0; i < myParams->n; i++) {
            myParams->t[i] = NTL::RandomBnd(myParams->q);
        }

        NTL::mat_ZZ td;
        _generateTrapdoor(myParams->A, td, myParams->n, myParams->l, myParams->q, myParams->beta);

        for (long i = 0; i < myParams->firstInputSize; i++) {
            for (long j = 0; j < myParams->firstInputSize; j++) {
                if (i > 2 || j > 2) {
                    continue;
                    myParams->u[i][j] = myParams->u[i % 3][j % 3];
                    continue;
                }
                _preSample(myParams->u[i][j], td, myParams->A,
                           myParams->v1[i] * v1Inv[j] * myParams->t, myParams->beta);
            }
            for (long j = 0; j < myParams->secondInputSize; j++) {
                if (i > 2 || j > 2) {
                    continue;
                    myParams->u[i][j + myParams->firstInputSize] =
                            myParams->u[i % 3][(j % 3) + myParams->firstInputSize];
                    continue;
                }
                _preSample(myParams->u[i][j + myParams->firstInputSize], td, myParams->A,
                           myParams->v1[i] * v2Inv[j] * myParams->t, myParams->beta);
            }
        }

        for (long i = 0; i < myParams->secondInputSize; i++) {
            for (long j = 0; j < myParams->firstInputSize; j++) {
                if (i > 2 || j > 2) {
                    continue;
                    myParams->u[i + myParams->firstInputSize][j] =
                            myParams->u[(i % 3) + myParams->firstInputSize][j % 3];
                    continue;
                }
                _preSample(myParams->u[i + myParams->firstInputSize][j], td, myParams->A,
                           myParams->v2[i] * v1Inv[j] * myParams->t, myParams->beta);
            }
            for (long j = 0; j < myParams->secondInputSize; j++) {
                if (i > 2 || j > 2) {
                    continue;
                    myParams->u[i + myParams->firstInputSize][j + myParams->firstInputSize] =
                            myParams->u[(i % 3) + myParams->firstInputSize][(j % 3) + myParams->firstInputSize];
                    continue;
                }
                _preSample(myParams->u[i + myParams->firstInputSize][j + myParams->firstInputSize], td, myParams->A,
                           myParams->v2[i] * v2Inv[j] * myParams->t, myParams->beta);
            }
        }

        auto end = chrono::steady_clock::now();
        auto ticks = chrono::duration_cast<chrono::milliseconds>(end - start).count();
        cout << "VC Setup ends. time: " << ticks << "ms" << endl;
    }

    void VectorCommitmentType2::commit(const MyFramework::VC::Commitment *commitment,
                                       const MyFramework::VC::Auxiliary *auxiliary,
                                       const MyFramework::VC::Params *params,
                                       const NTL::vec_ZZ &firstInput,
                                       const NTL::vec_ZZ &secondInput) {
        auto start = chrono::steady_clock::now();
        cout << "VC Commit starts" << endl;

        const auto myCommitment = (MyCommitment *) commitment;
        const auto myAuxiliary = (MyAuxiliary *) auxiliary;
        const auto myParams = (MyParams *) params;

        clear(myCommitment->c);
        myAuxiliary->x1 = firstInput;
        myAuxiliary->u1.SetLength(myParams->firstInputSize);
        for (long i = 0; i < myParams->firstInputSize; i++) {
            myCommitment->c += (myParams->isComplement ? NTL::InvMod(myParams->v1[i], myParams->q) : myParams->v1[i]) *
                    firstInput[i];
            myCommitment->c %= myParams->q;

            clear(myAuxiliary->u1[i]);
            myAuxiliary->u1[i].SetLength(myParams->l);
            for (long j = 0; j < myParams->firstInputSize; j++) {
                if (j == i) continue;
                myAuxiliary->u1[i] += firstInput[j] * (myParams->isComplement ? myParams->u[i % 3][j % 3] : myParams->u[j%3][i%3]);
            }
            for (long j = 0; j < myParams->secondInputSize; j++) {
                myAuxiliary->u1[i] += secondInput[j] * (myParams->isComplement
                                                            ? myParams->u[i%3][(j%3) + myParams->firstInputSize]
                                                            : myParams->u[(j%3) + myParams->firstInputSize][i%3]);
            }
        }

        myAuxiliary->x2 = secondInput;
        myAuxiliary->u2.SetLength(myParams->secondInputSize);
        for (long i = 0; i < myParams->secondInputSize; i++) {
            myCommitment->c += (myParams->isComplement ? NTL::InvMod(myParams->v2[i], myParams->q) : myParams->v2[i]) *
                    secondInput[i];
            myCommitment->c %= myParams->q;

            clear(myAuxiliary->u2[i]);
            myAuxiliary->u2[i].SetLength(myParams->l);
            for (long j = 0; j < myParams->firstInputSize; j++) {
                myAuxiliary->u2[i] += firstInput[j] * (myParams->isComplement
                                                           ? myParams->u[(i%3) + myParams->firstInputSize][j%3]
                                                           : myParams->u[j%3][(i%3) + myParams->firstInputSize]);
            }
            for (long j = 0; j < myParams->secondInputSize; j++) {
                if (j == i) continue;
                myAuxiliary->u2[i] += secondInput[j] *
                (myParams->isComplement
                     ? myParams->u[(i%3) + myParams->firstInputSize][(j%3) + myParams->firstInputSize]
                     : myParams->u[(j%3) + myParams->firstInputSize][(i%3) + myParams->firstInputSize]);
            }
        }

        auto end = chrono::steady_clock::now();
        auto ticks = chrono::duration_cast<chrono::milliseconds>(end - start).count();
        cout << "VC Commit ends. time: " << ticks << "ms" << endl;
    }

    void VectorCommitmentType2::open(const MyFramework::VC::OpeningProof *proof, const MyFramework::VC::Params *params,
                                     const MyFramework::VC::Auxiliary *auxiliary, const NTL::mat_ZZ &openingFunction1,
                                     const NTL::mat_ZZ &openingFunction2) {
        auto start = chrono::steady_clock::now();
        cout << "VC Open starts" << endl;

        const auto myProof = (MyOpeningProof *) proof;
        const auto myParams = (MyParams *) params;
        const auto myAuxiliary = (MyAuxiliary *) auxiliary;
        MyFramework::VC::Commitment *_commitmentTmp = new MyCommitment();
        MyFramework::VC::Auxiliary *_auxiliaryTmp = new MyAuxiliary();
        myProof->pi.SetLength(myParams->l);
        myProof->pi_eq.SetLength(myParams->l);
        myProof->pi_ip.SetLength(myParams->l);

        myProof->output = openingFunction1 * myAuxiliary->x1 + openingFunction2 * myAuxiliary->x2;

        NTL::vec_ZZ x1, x2;
        x1.SetLength(myParams->firstInputSize);
        x2.SetLength(myParams->secondInputSize);
        for (long i = 0; i < myParams->secondInputSize; i++) {
            x2[i] = myAuxiliary->x2[i] * myParams->h2[i];
        }
        myParams->isComplement = true;
        commit(_commitmentTmp, _auxiliaryTmp, myParams, x1, x2);
        const auto _commitment = (MyCommitment *) _commitmentTmp;
        const auto _auxiliary = (MyAuxiliary *) _auxiliaryTmp;
        myProof->_c = _commitment->c;
        myParams->isComplement = false;

        for (long i = 0; i < myParams->outputSize; i++) {
            for (long j = 0; j < myParams->firstInputSize; j++) {
                myProof->pi += myParams->h[i] * openingFunction1[i][j] * myAuxiliary->u1[j];
            }
            for (long j = 0; j < myParams->secondInputSize; j++) {
                myProof->pi += myParams->h[i] * openingFunction2[i][j] * myAuxiliary->u2[j];
            }
        }

        for (long i = 0; i < myParams->secondInputSize; i++) {
            myProof->pi_eq += myParams->h2[i] * NTL::InvMod(myParams->v2[i], myParams->q) * _auxiliary->u2[i];
        }

        for (long i = 0; i < myParams->secondInputSize; i++) {
            for (long j = 0; j < myParams->secondInputSize; j++) {
                if (i == j) continue;
                myProof->pi_ip += myAuxiliary->x2[i] * (x2[j] - myParams->h2[j]) * myParams->u[(i%3) + myParams->firstInputSize][(j%3) + myParams->firstInputSize];
            }
        }

        auto end = chrono::steady_clock::now();
        auto ticks = chrono::duration_cast<chrono::milliseconds>(end - start).count();
        cout << "VC Open ends. time: " << ticks << "ms" << endl;
    }

    bool VectorCommitmentType2::verify(const MyFramework::VC::Params *params, const NTL::mat_ZZ &openingFunction1,
                                       const NTL::mat_ZZ &openingFunction2,
                                       const MyFramework::VC::Commitment *commitment,
                                       const MyFramework::VC::OpeningProof *proof) {
        auto start = chrono::steady_clock::now();
        cout << "VC verify starts" << endl;

        const auto myParams = (MyParams *) params;
        const auto myCommitment = (MyCommitment *) commitment;
        const auto myProof = (MyOpeningProof *) proof;

        if (myProof->pi * myProof->pi > myParams->alpha) {
            return false;
        }

        if (myProof->pi_eq * myProof->pi_eq > myParams->alpha) {
            return false;
        }

        if (myProof->pi_ip * myProof->pi_ip > myParams->alpha) {
            return false;
        }

        auto end = chrono::steady_clock::now();
        auto ticks = chrono::duration_cast<chrono::milliseconds>(end - start).count();
        cout << "VC verify ends. time: " << ticks << "ms" << endl;

        return true;
    }
}
