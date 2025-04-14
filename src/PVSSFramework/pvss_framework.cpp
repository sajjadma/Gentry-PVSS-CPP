/* pvss_framework.cpp - a general PVSS framework
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
#include <pvss_framework.hpp>
#include <NTL/xdouble.h>
#include <NTL/ZZ_pX.h>

using namespace std;

namespace MyFramework {
    template<typename EncryptionType, typename VectorCommitmentType>
    void PVSS<EncryptionType, VectorCommitmentType>::setup(Params &params, const long securityParameter,
                                                           const long numberOfParties,
                                                           const long threshold) {
        params.numberOfParties = numberOfParties;
        params.threshold = threshold;

        NextPrime(params.prime, NTL::ZZ(numberOfParties + securityParameter));

        NTL::ZZ p3 = params.prime * params.prime * params.prime;

        this->encryptionSystem->setup(params.encryptionParams, securityParameter, p3);

        this->vectorCommitmentSystem->setup(params.vcParams, securityParameter,
                                            numberOfParties * params.encryptionParams->plainSize + threshold + 1,
                                            numberOfParties * params.encryptionParams->randomSize,
                                            numberOfParties * (params.encryptionParams->cipherSize + 1),
                                            p3,
                                            params.encryptionParams->randomBound,
                                            params.encryptionParams->coefficientBound > params.prime
                                                ? params.encryptionParams->coefficientBound
                                                : params.prime);
    }

    template<typename EncryptionType, typename VectorCommitmentType>
    void PVSS<EncryptionType, VectorCommitmentType>::generateKey(Encryption::KeyPair &key, const Params &params) {
        this->encryptionSystem->generateKey(&key, params.encryptionParams);
    }

    template<typename EncryptionType, typename VectorCommitmentType>
    bool PVSS<EncryptionType, VectorCommitmentType>::verifyKey(const Params &params,
                                                               const Encryption::PublicKey *publicKey,
                                                               const Encryption::KeyProof *proof) {
        return this->encryptionSystem->verifyKey(params.encryptionParams, publicKey, proof);
    }

    template<typename EncryptionType, typename VectorCommitmentType>
    void PVSS<EncryptionType, VectorCommitmentType>::distribute(DistributionProof &proof, const Params &params,
                                                                const vector<Encryption::PublicKey *> &publicKeys,
                                                                const NTL::ZZ &secret) {
        NTL::vec_ZZ firstInput, secondInput;
        firstInput.SetLength(params.vcParams->firstInputSize);
        secondInput.SetLength(params.vcParams->secondInputSize);
        NTL::mat_ZZ M1, M2;
        M1.SetDims(params.vcParams->outputSize, params.vcParams->firstInputSize);
        M2.SetDims(params.vcParams->outputSize, params.vcParams->secondInputSize);
        proof.encryptedShares.clear();

        // Polynomial part starts
        firstInput[0] = secret;
        for (long i = 1; i <= params.threshold; i++) {
            RandomBnd(firstInput[i], params.prime);
        }
        // Polynomial part ends

        long outputIndex = 0;
        for (long i = 0; i < params.numberOfParties; i++) {
            M1.put(outputIndex, 0, NTL::to_ZZ(1));
            NTL::ZZ share = firstInput[0];
            for (long j = 1; j < params.threshold + 1; j++) {
                M1.put(outputIndex, j, ((i + 1) * M1[outputIndex][j - 1]) % params.prime);
                share = share + firstInput[j] * M1[outputIndex][j];
            }

            NTL::vec_ZZ encodedShare;
            encodedShare.SetLength(params.encryptionParams->plainSize);
            for (long j = 0; j < params.encryptionParams->plainSize; j++) {
                const long index = params.threshold + 1 + (i * params.encryptionParams->plainSize) + j;
                encodedShare[j] = share % params.encryptionParams->plainBound;
                firstInput[index] = encodedShare[j];
                M1.put(outputIndex, index,
                       j == 0 ? NTL::to_ZZ(1) : params.encryptionParams->plainBound * M1[outputIndex][index - 1]);
                share = share / params.encryptionParams->plainBound;
            }

            NTL::vec_ZZ r;
            NTL::ZZ tmpR;
            r.SetLength(params.encryptionParams->randomSize);
            for (long j = 0; j < params.encryptionParams->randomSize; j++) {
                RandomBnd(tmpR, params.encryptionParams->randomBound);
                r[j] = tmpR;
                for (long k = 0; k < params.vcParams->secondInputPartitionSize; k++) {
                    const long index = (i * params.encryptionParams->randomSize) +
                                       (j * params.vcParams->secondInputPartitionSize) + k;
                    secondInput[index] = tmpR % params.vcParams->secondInputBound;
                    tmpR = tmpR / params.vcParams->secondInputBound;
                }
            }

            NTL::mat_ZZ f1, f2;
            this->encryptionSystem->generateEncryptionFunctionFromInput(f1, f2, params.encryptionParams, publicKeys[i],
                                                                        encodedShare, r);

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
                                     M1[outputIndex][index - 1]);
                    }
                }
            }
        }

        VC::Auxiliary *auxiliary;
        this->vectorCommitmentSystem->commit(proof.commitment, auxiliary, params.vcParams, firstInput, secondInput);
        this->vectorCommitmentSystem->open(proof.proof, params.vcParams, auxiliary, M1, M2);
    }

    template<typename EncryptionType, typename VectorCommitmentType>
    bool PVSS<EncryptionType, VectorCommitmentType>::verifyDistribution(
        const Params &params, const vector<Encryption::PublicKey *> &publicKeys,
        const DistributionProof &proof) {
        NTL::mat_ZZ M1, M2;
        M1.SetDims(params.vcParams->outputSize, params.vcParams->firstInputSize);
        M2.SetDims(params.vcParams->outputSize, params.vcParams->secondInputSize);

        long outputIndex = 0;
        for (long i = 0; i < params.numberOfParties; i++) {
            M1.put(outputIndex, 0, NTL::to_ZZ(1));
            for (long j = 1; j < params.threshold + 1; j++) {
                if (proof.proof->output[outputIndex] != 0) {
                    return false;
                }

                M1.put(outputIndex, j, ((i + 1) * M1[outputIndex][j - 1]) % params.prime);
            }

            for (long j = 0; j < params.encryptionParams->plainSize; j++) {
                const long index = params.threshold + 1 + (i * params.encryptionParams->plainSize) + j;
                M1.put(outputIndex, index,
                       j == 0 ? NTL::to_ZZ(1) : params.encryptionParams->plainBound * M1[outputIndex][index - 1]);
            }

            NTL::mat_ZZ f1, f2;
            this->encryptionSystem->generateEncryptionFunctionFromOutput(
                f1, f2, params.encryptionParams, publicKeys[i],
                proof.encryptedShares[i]);

            for (long j = 0; j < params.encryptionParams->cipherSize; j++) {
                outputIndex++;
                if (proof.proof->output[outputIndex] != proof.encryptedShares[i][j]) {
                    return false;
                }

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
                                     M1[outputIndex][index - 1]);
                    }
                }
            }
        }

        return this->vectorCommitmentSystem->verify(params.vcParams, M1, M2, proof.commitment, proof.proof);
    }

    template<typename EncryptionType, typename VectorCommitmentType>
    void PVSS<EncryptionType, VectorCommitmentType>::decryptShare(DecryptionProof &proof, const Params &params,
                                                                  const Encryption::PublicKey *publicKey,
                                                                  const Encryption::PrivateKey *privateKey,
                                                                  const NTL::vec_ZZ &encryptedShare) {
        this->encryptionSystem->decrypt(proof.proof, params.encryptionParams, publicKey, privateKey, encryptedShare);
        proof.decryptedShare = 0;
        NTL::ZZ pow = NTL::to_ZZ(1);
        for (long i = 0; i < params.encryptionParams->plainSize; i++) {
            proof.decryptedShare += proof.proof->decryptedValues[i] * pow;
            pow *= params.encryptionParams->plainBound;
        }
    }

    template<typename EncryptionType, typename VectorCommitmentType>
    bool PVSS<EncryptionType, VectorCommitmentType>::verifyDecryption(const Params &params,
                                                                      const Encryption::PublicKey *publicKey,
                                                                      const NTL::vec_ZZ &encryptedShare,
                                                                      const DecryptionProof &proof) {
        NTL::ZZ share = NTL::ZZ::zero();
        NTL::ZZ pow = NTL::to_ZZ(1);
        for (long i = 0; i < params.encryptionParams->plainSize; i++) {
            share += proof.proof->decryptedValues[i] * pow;
            pow *= params.encryptionParams->plainBound;
        }
        return share == proof.decryptedShare && this->encryptionSystem->verifyDecryption(
                   params.encryptionParams, publicKey, encryptedShare, proof.proof);
    }

    // TODO: Replace with lagrange interpolation to improve performance
    template<typename EncryptionType, typename VectorCommitmentType>
    void PVSS<EncryptionType, VectorCommitmentType>::reconstruct(NTL::ZZ &reconstruction, const Params &params,
                                                                 const vector<NTL::ZZ> &decryptedShares) {
        NTL::ZZ_pPush push(params.prime);
        NTL::vec_ZZ_p a, b;
        a.SetLength(params.numberOfParties);
        b.SetLength(params.numberOfParties);
        for (long i = 0; i < params.numberOfParties; i++) {
            a[i] = NTL::to_ZZ_p(i + 1);
            b[i] = NTL::to_ZZ_p(decryptedShares[i]);
        }

        const NTL::ZZ_pX f = NTL::interpolate(a, b);
        reconstruction = NTL::rep(NTL::eval(f, NTL::ZZ_p::zero()));
    }

    template class PVSS<MyEncryption::EncryptionType1, MyVectorCommitment::VectorCommitmentType1>;
}
