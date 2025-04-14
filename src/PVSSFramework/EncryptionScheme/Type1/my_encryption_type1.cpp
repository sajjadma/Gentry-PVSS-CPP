/* my_encryption_type1.cpp
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
#include <NTL/xdouble.h>
#include <sodium/crypto_generichash.h>

using namespace std;

namespace MyEncryption {
    void EncryptionType1::_generateTrapdoor(NTL::mat_ZZ &A, NTL::mat_ZZ &trapdoor, const long n, const long m,
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

    void EncryptionType1::_preSample(NTL::vec_ZZ &x, const NTL::mat_ZZ &trapdoor, const NTL::mat_ZZ &A,
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

    void EncryptionType1::_hash(NTL::vec_ZZ &hash, const NTL::mat_ZZ &A) {
        hash.SetLength(A.NumRows());
        unsigned char _hash[A.NumRows()];

        crypto_generichash_state state;
        crypto_generichash_init(&state, (unsigned char *) "Key Hash", 8, sizeof(hash));
        for (long i = 0; i < A.NumRows(); i++) {
            for (long j = 0; j < A.NumCols(); j++) {
                const long size = NumBytes(A[i][j]);
                unsigned char buf[size];
                NTL::BytesFromZZ(buf, A[i][j], size);
                crypto_generichash_update(&state, buf, size);
            }
        }
        crypto_generichash_final(&state, _hash, sizeof(_hash));

        for (long i = 0; i < A.NumRows(); i++) {
            hash[i] = _hash[i];
        }
    }

    // TODO: Set Parameters Properly
    void EncryptionType1::setup(const MyFramework::Encryption::Params *params, const long securityParameter,
                                const NTL::ZZ plainBound) {
        const auto myParams = (MyParams *) params;
        const long k = NTL::NumBits(plainBound);
        myParams->module = NTL::power2_ZZ(k);
        myParams->l = 8; // TODO: یا هر مقدار مناسب دیگر
        myParams->m = 2 * myParams->l * k;
        myParams->d = k; // TODO: یا هر مقدار مناسب دیگر
        myParams->randomBound = NTL::power2_ZZ(k / 2 - 2);
        myParams->plainBound = myParams->module;
        myParams->coefficientBound = myParams->module;
        myParams->plainSize = 1;
        myParams->randomSize = myParams->m * myParams->d + myParams->l + myParams->m + myParams->d;
        myParams->cipherSize = myParams->l * myParams->d + myParams->m + myParams->d;
    }

    void EncryptionType1::generateKey(const MyFramework::Encryption::KeyPair *key,
                                      const MyFramework::Encryption::Params *params) {
        const auto myParams = (MyParams *) params;
        const auto privateKey = (MyPrivateKey *) key->privateKey;
        const auto publicKey = (MyPublicKey *) key->publicKey;
        const auto proof = (MyKeyProof *) key->proof;

        NTL::vec_ZZ b;
        _generateTrapdoor(publicKey->A, privateKey->trapdoor, myParams->l, myParams->m, myParams->module,
                          myParams->randomBound);
        _hash(b, publicKey->A);
        _preSample(proof->x, privateKey->trapdoor, publicKey->A, b, myParams->randomBound);
    }

    bool EncryptionType1::verifyKey(const MyFramework::Encryption::Params *params,
                                    const MyFramework::Encryption::PublicKey *publicKey,
                                    const MyFramework::Encryption::KeyProof *proof) {
        const auto myParams = (MyParams *) params;
        const auto myPublicKey = (MyPublicKey *) publicKey;
        const auto myProof = (MyKeyProof *) proof;
        NTL::vec_ZZ b, y = myPublicKey->A * myProof->x;
        _hash(b, myPublicKey->A);

        for (long i = 0; i < b.length(); ++i) {
            if (b[i] % myParams->module != y[i] % myParams->module) {
                return false;
            }
        }

        return true;
    }

    void EncryptionType1::generateEncryptionFunctionFromInput(NTL::mat_ZZ &f1, NTL::mat_ZZ &f2,
                                                              const MyFramework::Encryption::Params *params,
                                                              const MyFramework::Encryption::PublicKey *publicKey,
                                                              const NTL::vec_ZZ &plainValues,
                                                              const NTL::vec_ZZ &randomValues) {
        const auto myParams = (MyParams *) params;
        const auto myPublicKey = (MyPublicKey *) publicKey;
        NTL::clear(f1);
        f1.SetDims(myParams->cipherSize, 1);
        NTL::clear(f2);
        f2.SetDims(myParams->cipherSize, myParams->randomSize);
        NTL::mat_ZZ B;
        B.SetDims(myParams->m, myParams->d);
        long randomIndex = 0;
        long outputIndex = 0;

        for (long i = myParams->cipherSize - myParams->d; i < myParams->cipherSize; i++) {
            if (i == myParams->cipherSize - myParams->d) {
                f1[i][0] = 1;
            } else {
                f1[i][0] = f1[i - 1][0] * 2; // TODO: If param.d changes 2 must change too
            }
        }

        for (long j = 0; j < myParams->d; j++) {
            for (long i = 0; i < myParams->m; i++) {
                B[i][j] = randomValues[randomIndex];
                randomIndex++;
            }
        }

        outputIndex = myParams->l * myParams->d;
        NTL::mat_ZZ U = myPublicKey->A * B;
        for (long i = 0; i < myParams->l; i++) {
            for (long j = 0; j < myParams->m; j++) {
                for (long n = 0; n < myParams->d; n++) {
                    // f2 coefficient corresponding with B
                    f2[n * myParams->l + i][n * myParams->m + j] = myPublicKey->A[i][j];
                }

                // f2 coefficient corresponding with f
                f2[outputIndex + j][randomIndex + i] = myPublicKey->A[i][j];
            }
        }

        for (long i = 0; i < myParams->m; i++) {
            // f2 coefficient corresponding with e
            f2[outputIndex + i][randomIndex + myParams->l + i] = 1;
        }

        outputIndex += myParams->m;
        for (long i = 0; i < myParams->l; i++) {
            for (long j = 0; j < myParams->d; j++) {
                // f2 coefficient corresponding with f
                f2[outputIndex + j][randomIndex + i] = U[i][j];
            }
        }

        for (long i = 0; i < myParams->d; i++) {
            // f2 coefficient corresponding with e^t
            f2[outputIndex + i][randomIndex + myParams->l + myParams->m + i] = 1;
        }
    }

    void EncryptionType1::generateEncryptionFunctionFromOutput(NTL::mat_ZZ &f1, NTL::mat_ZZ &f2,
                                                               const MyFramework::Encryption::Params *params,
                                                               const MyFramework::Encryption::PublicKey *publicKey,
                                                               const NTL::vec_ZZ &cipherValues) {
        const auto myParams = (MyParams *) params;
        const auto myPublicKey = (MyPublicKey *) publicKey;
        NTL::clear(f1);
        f1.SetDims(myParams->cipherSize, 1);
        NTL::clear(f2);
        f2.SetDims(myParams->cipherSize, myParams->randomSize);
        NTL::mat_ZZ U;
        U.SetDims(myParams->l, myParams->d);
        long randomIndex = 0;
        long outputIndex = 0;

        for (long i = myParams->cipherSize - myParams->d; i < myParams->cipherSize; i++) {
            if (i == myParams->cipherSize - myParams->d) {
                f1[i][0] = 1;
            } else {
                f1[i][0] = f1[i - 1][0] * 2; // TODO: If param.d changes 2 must change too
            }
        }

        for (long j = 0; j < myParams->d; j++) {
            for (long i = 0; i < myParams->l; i++) {
                U[i][j] = cipherValues[outputIndex];
                outputIndex++;
            }
        }

        randomIndex = myParams->m * myParams->d;
        for (long i = 0; i < myParams->l; i++) {
            for (long j = 0; j < myParams->m; j++) {
                for (long n = 0; n < myParams->d; n++) {
                    // f2 coefficient corresponding with B
                    f2[n * myParams->l + i][n * myParams->m + j] = myPublicKey->A[i][j];
                }

                // f2 coefficient corresponding with f
                f2[outputIndex + j][randomIndex + i] = myPublicKey->A[i][j];
            }
        }

        for (long i = 0; i < myParams->m; i++) {
            // f2 coefficient corresponding with e
            f2[outputIndex + i][randomIndex + myParams->l + i] = 1;
        }

        outputIndex += myParams->m;
        for (long i = 0; i < myParams->l; i++) {
            for (long j = 0; j < myParams->d; j++) {
                // f2 coefficient corresponding with f
                f2[outputIndex + j][randomIndex + i] = U[i][j];
            }
        }

        for (long i = 0; i < myParams->d; i++) {
            // f2 coefficient corresponding with e^t
            f2[outputIndex + i][randomIndex + myParams->l + myParams->m + i] = 1;
        }
    }

    void EncryptionType1::decrypt(const MyFramework::Encryption::DecryptionProof *proof,
                                  const MyFramework::Encryption::Params *params,
                                  const MyFramework::Encryption::PublicKey *publicKey,
                                  const MyFramework::Encryption::PrivateKey *privateKey,
                                  const NTL::vec_ZZ &cipherValues) {
        const auto myProof = (MyDecryptionProof *) proof;
        const auto myParams = (MyParams *) params;
        const auto myPublicKey = (MyPublicKey *) publicKey;
        const auto myPrivateKey = (MyPrivateKey *) privateKey;
        NTL::vec_ZZ tmpU, tmpB, h, C;
        tmpU.SetLength(myParams->l);
        h.SetLength(myParams->m);
        C.SetLength(myParams->d);
        myProof->B.SetDims(myParams->m, myParams->d);
        myProof->decryptedValues.SetLength(1);
        myProof->decryptedValues[0] = 0;
        long outputIndex = 0;

        for (long j = 0; j < myParams->d; j++) {
            for (long i = 0; i < myParams->l; i++) {
                tmpU[i] = cipherValues[outputIndex];
                outputIndex++;
            }

            _preSample(tmpB, myPrivateKey->trapdoor, myPublicKey->A, tmpU, myParams->randomBound);
            for (long i = 0; i < myParams->m; i++) {
                myProof->B[i][j] = tmpB[i];
            }
        }

        for (long i = 0; i < myParams->m; i++) {
            h[i] = cipherValues[outputIndex];
            outputIndex++;
        }

        for (long i = 0; i < myParams->d; i++) {
            C[i] = cipherValues[outputIndex];
            outputIndex++;
        }

        NTL::vec_ZZ tmpM = C - (h * myProof->B);
        const NTL::ZZ bound = myParams->module / 2;
        // TODO: If param.d changes 2 must change too
        NTL::ZZ pow = myParams->module / 2;
        for (long i = myParams->d - 1; i >= 0; i--) {
            if (((tmpM[i] - pow * myProof->decryptedValues[0]) % myParams->module) > bound) {
                myProof->decryptedValues[0] += bound / pow;
            }

            tmpM[i] = pow * myProof->decryptedValues[0];
            pow /= 2;
        }

        myProof->e = C - h * myProof->B - tmpM;
    }

    bool EncryptionType1::verifyDecryption(const MyFramework::Encryption::Params *params,
                                           const MyFramework::Encryption::PublicKey *publicKey,
                                           const NTL::vec_ZZ &cipherValues,
                                           const MyFramework::Encryption::DecryptionProof *proof) {
        const auto myParams = (MyParams *) params;
        const auto myPublicKey = (MyPublicKey *) publicKey;
        const auto myProof = (MyDecryptionProof *) proof;

        if (myProof->e * myProof->e > myParams->module / 4) {
            return false;
        }

        for (long i = 0; i < myParams->m; i++) {
            if (myProof->B[i] * myProof->B[i] > myParams->module / 4) {
                return false;
            }
        }

        NTL::vec_ZZ h;
        NTL::mat_ZZ U = myPublicKey->A * myProof->B;
        h.SetLength(myParams->m);
        long outputIndex = 0;

        for (long j = 0; j < myParams->d; j++) {
            for (long i = 0; i < myParams->l; i++) {
                if (U[i][j] != cipherValues[outputIndex]) {
                    return false;
                }
                outputIndex++;
            }
        }

        for (long i = 0; i < myParams->m; i++) {
            h[i] = cipherValues[outputIndex];
            outputIndex++;
        }

        NTL::vec_ZZ C = h * myProof->B + myProof->e;
        NTL::ZZ pow;
        set(pow);
        for (long i = 0; i < myParams->d; i++) {
            if (C[i] + (pow * myProof->decryptedValues[0]) != cipherValues[outputIndex]) {
                return false;
            }
            outputIndex++;
            pow *= 2; // TODO: If param.d changes 2 must change too
        }

        return true;
    }
}
