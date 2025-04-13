/* my_vector_commitment_type1.cpp
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
#include <NTL/ZZ_pX.h>
#include <sodium/crypto_generichash.h>

using namespace std;

namespace MyVectorCommitment {
    void VectorCommitmentType1::_generateTrapdoor(NTL::mat_ZZ &A, NTL::mat_ZZ &trapdoor, const long n, const long m,
                                                  const NTL::ZZ &module, const NTL::ZZ &bound) {
        const long k = NTL::NumBits(module);
        if (m < n * k) throw std::invalid_argument("EncryptionType1::_generateTrapdoor");

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
                G[i][index] = j == 0 ? 1 : G[i][index - 1] * 2;
            }
        }

        G = G - A_hat * trapdoor;

        for (long i = 0; i < n; i++) {
            for (long j = 0; j < _n; j++) {
                A[i][_m + j] = G[i][j];
            }
        }
    }

    void VectorCommitmentType1::_preSample(NTL::vec_ZZ &x, const NTL::mat_ZZ &trapdoor, const NTL::mat_ZZ &A,
                                           const NTL::vec_ZZ &b, const NTL::ZZ &bound) {
        const long n = A.NumRows(), m = trapdoor.NumCols();
        const long k = m / n;
        if (b.length() != n || trapdoor.NumRows() != n || m != n * k)
            throw std::invalid_argument("EncryptionType1::_preSample");

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
    void VectorCommitmentType1::setup(MyFramework::VC::Params *params, long securityParameter, long firstInputSize,
                                      long secondInputSize, long outputSize, NTL::ZZ firstInputBound,
                                      NTL::ZZ secondInputBound, NTL::ZZ coefficientBound) {
        MyParams myParams;
        const long n = NTL::NumBits(secondInputBound);
        const long k = NTL::NumBits(
            (firstInputBound > coefficientBound ? firstInputBound : coefficientBound) * firstInputSize * outputSize);
        myParams.firstInputSize = firstInputSize;
        myParams.secondInputSize = secondInputSize * n;
        myParams.secondInputPartitionSize = n;
        myParams.outputSize = outputSize;
        myParams.firstInputBound = NTL::power2_ZZ(k);
        myParams.secondInputBound = 2;
        myParams.coefficientBound = myParams.firstInputBound;
        NTL::vec_ZZ v1Inv, v2Inv;
        v1Inv.SetLength(myParams.firstInputSize);
        myParams.v1.SetLength(myParams.firstInputSize);
        v2Inv.SetLength(myParams.secondInputSize);
        myParams.v2.SetLength(myParams.secondInputSize);
        myParams.h.SetLength(myParams.outputSize);
        myParams.h2.SetLength(myParams.secondInputSize);
        myParams.u.SetDims(myParams.firstInputSize + myParams.secondInputSize,
                           myParams.firstInputSize + myParams.secondInputSize);
        myParams.n = 16 * k;
        myParams.l = 32 * myParams.n * k;
        myParams.p = NTL::power2_ZZ(4 * k);
        myParams.q = NTL::power2_ZZ(16 * k);
        myParams.beta = NTL::power2_ZZ(k);
        myParams.alpha = NTL::power2_ZZ(12 * k);
        myParams.t.SetLength(myParams.n);

        for (long i = 0; i < myParams.firstInputSize; i++) {
            while (true) {
                try {
                    v1Inv[i] = NTL::RandomBnd(myParams.q);
                    NTL::InvMod(myParams.v1[i], v1Inv[i], myParams.q);
                    break;
                } catch ([[maybe_unused]] exception &e) {
                }
            }
        }

        for (long i = 0; i < myParams.secondInputSize; i++) {
            while (true) {
                try {
                    v2Inv[i] = NTL::RandomBnd(myParams.p);
                    NTL::InvMod(myParams.v2[i], v2Inv[i], myParams.q);
                    break;
                } catch ([[maybe_unused]] exception &e) {
                }
            }
        }

        for (long i = 0; i < myParams.outputSize; i++) {
            myParams.h[i] = NTL::RandomBnd(myParams.p);
        }

        for (long i = 0; i < myParams.secondInputSize; i++) {
            myParams.h2[i] = NTL::RandomBnd(myParams.p);
        }

        for (long i = 0; i < myParams.n; i++) {
            myParams.t[i] = NTL::RandomBnd(myParams.q);
        }

        NTL::mat_ZZ td;
        _generateTrapdoor(myParams.A, td, myParams.n, myParams.l, myParams.q, myParams.beta);

        for (long i = 0; i < myParams.firstInputSize; i++) {
            for (long j = 0; j < myParams.firstInputSize; j++) {
                if (i == j) continue;
                _preSample(myParams.u[i][j], td, myParams.A,
                           myParams.v1[i] * v1Inv[j] * myParams.t, myParams.beta);
            }
            for (long j = 0; j < myParams.secondInputSize; j++) {
                _preSample(myParams.u[i][j + myParams.firstInputSize], td, myParams.A,
                           myParams.v1[i] * v2Inv[j] * myParams.t, myParams.beta);
            }
        }

        for (long i = 0; i < myParams.secondInputSize; i++) {
            for (long j = 0; j < myParams.firstInputSize; j++) {
                _preSample(myParams.u[i + myParams.firstInputSize][j], td, myParams.A,
                           myParams.v2[i] * v1Inv[j] * myParams.t, myParams.beta);
            }
            for (long j = 0; j < myParams.secondInputSize; j++) {
                if (i == j) continue;
                _preSample(myParams.u[i + myParams.firstInputSize][j + myParams.firstInputSize], td, myParams.A,
                           myParams.v2[i] * v2Inv[j] * myParams.t, myParams.beta);
            }
        }

        params = &myParams;
    }

    void VectorCommitmentType1::commit(MyFramework::VC::Commitment *commitment, MyFramework::VC::Auxiliary *auxiliary,
                                       const MyFramework::VC::Params *params, const NTL::vec_ZZ &firstInput,
                                       const NTL::vec_ZZ &secondInput) {
        MyCommitment myCommitment;
        MyAuxiliary myAuxiliary;
        const auto myParams = (MyParams *) params;

        clear(myCommitment.c);
        myAuxiliary.x1 = firstInput;
        myAuxiliary.u1.SetLength(myParams->firstInputSize);
        for (long i = 0; i < myParams->firstInputSize; i++) {
            myCommitment.c += (myParams->isComplement ? NTL::InvMod(myParams->v1[i], myParams->q) : myParams->v1[i]) *
                    firstInput[i];
            myCommitment.c %= myParams->q;

            clear(myAuxiliary.u1[i]);
            myAuxiliary.u1[i].SetLength(myParams->l);
            for (long j = 0; j < myParams->firstInputSize; j++) {
                if (j == i) continue;
                myAuxiliary.u1[i] += firstInput[j] * (myParams->isComplement ? myParams->u[i][j] : myParams->u[j][i]);
            }
            for (long j = 0; j < myParams->secondInputSize; j++) {
                myAuxiliary.u1[i] += secondInput[j] * (myParams->isComplement
                                                           ? myParams->u[i][j + myParams->firstInputSize]
                                                           : myParams->u[j + myParams->firstInputSize][i]);
            }
        }

        myAuxiliary.x2 = secondInput;
        myAuxiliary.u2.SetLength(myParams->secondInputSize);
        for (long i = 0; i < myParams->secondInputSize; i++) {
            myCommitment.c += (myParams->isComplement ? NTL::InvMod(myParams->v2[i], myParams->q) : myParams->v2[i]) *
                    secondInput[i];
            myCommitment.c %= myParams->q;

            clear(myAuxiliary.u2[i]);
            myAuxiliary.u2[i].SetLength(myParams->l);
            for (long j = 0; j < myParams->firstInputSize; j++) {
                myAuxiliary.u2[i] += firstInput[j] * (myParams->isComplement
                                                          ? myParams->u[i + myParams->firstInputSize][j]
                                                          : myParams->u[j][i + myParams->firstInputSize]);
            }
            for (long j = 0; j < myParams->secondInputSize; j++) {
                if (j == i) continue;
                myAuxiliary.u2[i] += secondInput[j] *
                (myParams->isComplement
                     ? myParams->u[i + myParams->firstInputSize][j + myParams->firstInputSize]
                     : myParams->u[j + myParams->firstInputSize][i + myParams->firstInputSize]);
            }
        }

        commitment = &myCommitment;
        auxiliary = &myAuxiliary;
    }

    void VectorCommitmentType1::open(MyFramework::VC::OpeningProof *proof, const MyFramework::VC::Params *params,
                                     MyFramework::VC::Auxiliary *auxiliary, const NTL::mat_ZZ &openingFunction1,
                                     const NTL::mat_ZZ &openingFunction2) {
        MyOpeningProof myProof;
        auto myParams = (MyParams *) params;
        auto myAuxiliary = (MyAuxiliary *) auxiliary;
        MyCommitment _commitment;
        MyAuxiliary _auxiliary;
        myProof.pi.SetLength(myParams->l);
        myProof.pi_eq.SetLength(myParams->l);
        myProof.pi_ip.SetLength(myParams->l);

        myProof.output = openingFunction1 * myAuxiliary->x1 + openingFunction2 * myAuxiliary->x2;

        NTL::vec_ZZ x1, x2;
        x1.SetLength(myParams->firstInputSize);
        x2.SetLength(myParams->secondInputSize);
        for (long i = 0; i < myParams->secondInputSize; i++) {
            x2[i] = myAuxiliary->x2[i] * myParams->h2[i];
        }
        myParams->isComplement = true;
        commit(&_commitment, &_auxiliary, myParams, x1, x2);
        myProof._c = _commitment.c;
        myParams->isComplement = false;

        for (long i = 0; i < myParams->outputSize; i++) {
            for (long j = 0; j < myParams->firstInputSize; j++) {
                myProof.pi += myParams->h[i] * openingFunction1[i][j] * myAuxiliary->u1[j];
            }
            for (long j = 0; j < myParams->secondInputSize; j++) {
                myProof.pi += myParams->h[i] * openingFunction2[i][j] * myAuxiliary->u2[j];
            }
        }

        for (long i = 0; i < myParams->secondInputSize; i++) {
            myProof.pi_eq += myParams->h2[i] * NTL::InvMod(myParams->v2[i], myParams->q) * _auxiliary.u2[i];
        }

        for (long i = 0; i < myParams->secondInputSize; i++) {
            for (long j = 0; j < myParams->secondInputSize; j++) {
                if (i == j) continue;
                myProof.pi_ip += myAuxiliary->x2[i] * (x2[j] - myParams->h2[j]) * myParams->u[i][j];
            }
        }

        proof = &myProof;
    }

    bool VectorCommitmentType1::verify(const MyFramework::VC::Params *params, const NTL::mat_ZZ &openingFunction1,
                                       const NTL::mat_ZZ &openingFunction2,
                                       const MyFramework::VC::Commitment *commitment,
                                       const MyFramework::VC::OpeningProof *proof) {
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

        return true;
    }
}
