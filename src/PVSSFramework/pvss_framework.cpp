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

#include <pvss_framework.hpp>

using namespace std;

namespace MyFramework {
    template<typename R>
    Params PVSS<R>::setup(int securityParameter, int numberOfParties, int threshold) {
        auto params = Params();
        params.numberOfParties = numberOfParties;
        params.threshold = threshold;

        params.prime = numberOfParties + 1; // TODO: select prime

        params.encryptionParams = this->encryptionSystem->setup(securityParameter);

        int inputSize = numberOfParties * (params.encryptionParams.randomSize + 1) + threshold + 1;
        int outputSize = numberOfParties * (params.encryptionParams.cipherTextSize + 1);
        params.vcParams = this->vectorCommitmentSystem->setup(securityParameter, inputSize, outputSize,
                                                              params.encryptionParams.randomSize,
                                                              params.encryptionParams.randomBound);

        return params;
    }

    template<typename R>
    Encryption::KeyPair PVSS<R>::generateKey(const Params &params, int index) {
        return this->encryptionSystem->generateKey(params.encryptionParams);
    }

    template<typename R>
    bool PVSS<R>::verifyKey(const Params &params, int index, const Encryption::PublicKey &publicKey,
                            const Encryption::KeyProof &proof) {
        return this->encryptionSystem->verifyKey(params.encryptionParams, publicKey, proof);
    }

    template<typename R>
    DistributionProof<R> PVSS<R>::distribute(const Params &params,
                                             const vector<Encryption::PublicKey> &publicKeys, const R &secret) {
        vector<R> polynomial(params.threshold + 1);
        polynomial[0] = secret;
        for (int i = 1; i <= params.threshold; i++) {
            polynomial[i] = params.prime; // TODO: random
        }

        vector<R> x(polynomial);
        vector<vector<R> > M;
        vector<vector<R> > encryptedShares(params.numberOfParties);
        for (int i = 0; i < params.numberOfParties; i++) {
            vector<R> b(params.threshold + 1);
            b[0] = 1;
            R share = polynomial[0];
            for (int j = 1; j < params.threshold + 1; j++) {
                b[j] = (i * b[j - 1]) % params.prime;
                share += polynomial[j] * b[j];
            }

            x.emplace_back(share);
            M.emplace_back(vector<R>(params.vcParams.inputSize));
            copy(b.begin(), b.end(), M.back().begin());
            M.back()[x.size() - 1] = -1;

            vector<R> r(params.encryptionParams.randomSize);
            for (int j = 0; j < params.encryptionParams.randomSize; j++) {
                r[j] = params.encryptionParams.randomBound; // TODO: random
            }

            vector<vector<R> > encryptionFunction = this->encryptionSystem->generateEncryptionFunctionFromInput(
                params.encryptionParams, publicKeys[i], share, r);
            vector<R> a(1);
            a[0] = share;
            a.insert(a.end(), r.begin(), r.end());
            encryptedShares[i] = encryptionFunction * a;

            x.insert(x.end(), r.begin(), r.end());
            for (int j = 0; j < params.encryptionParams.cipherTextSize; j++) {
                M.emplace_back(vector<R>(params.vcParams.inputSize));
                copy(encryptionFunction[j].begin(), encryptionFunction[j].end(),
                     M.back().begin() + x.size() - a.size());
            }
        }

        pair<VC::Commitment, VC::Auxiliary> commitments = this->vectorCommitmentSystem->commit(
            params.vcParams, x);
        VC::OpeningProof<R> openingProof = this->vectorCommitmentSystem->open(params.vcParams, commitments.second, M);

        DistributionProof<R> proof = DistributionProof<R>();
        proof.commitment = commitments.first;
        proof.proof = openingProof;
        proof.encryptedShares = encryptedShares;

        return proof;
    }

    template<typename R>
    bool PVSS<R>::verifyDistribution(const Params &params, const vector<Encryption::PublicKey> &publicKeys,
                                     const DistributionProof<R> &proof) {
        vector<vector<R> > M;
        int index = params.threshold + 1;
        for (int i = 0; i < params.numberOfParties; i++) {
            if (proof.proof.output[M.size()] != 0) {
                return false;
            }

            vector<R> b(params.threshold + 1);
            b[0] = 1;
            for (int j = 1; j < params.threshold + 1; j++) {
                b[j] = (i * b[j - 1]) % params.prime;
            }

            ++index;
            M.emplace_back(vector<R>(params.vcParams.inputSize));
            copy(b.begin(), b.end(), M.back().begin());
            M.back()[index - 1] = -1;

            vector<vector<R> > encryptionFunction = this->encryptionSystem->generateEncryptionFunctionFromOutput(
                params.encryptionParams, publicKeys[i], proof.encryptedShares[i]);

            index += params.encryptionParams.randomSize;
            for (int j = 0; j < params.encryptionParams.cipherTextSize; j++) {
                if (proof.proof.output[M.size()] != proof.encryptedShares[i][j]) {
                    return false;
                }

                M.emplace_back(vector<R>(params.vcParams.inputSize));
                copy(encryptionFunction[j].begin(), encryptionFunction[j].end(),
                     M.back().begin() + index - params.encryptionParams.randomSize - 1);
            }
        }

        return this->vectorCommitmentSystem->verify(params.vcParams, M, proof.commitment, proof.proof);
    }

    template<typename R>
    Encryption::DecryptionProof<R> PVSS<R>::decryptShare(const Params &params, int index,
                                                         const Encryption::PrivateKey &privateKey,
                                                         const vector<R> &encryptedShare) {
        return this->encryptionSystem->decrypt(params.encryptionParams, privateKey, encryptedShare);
    }

    template<typename R>
    bool PVSS<R>::verifyDecryption(const Params &params, int index, const Encryption::PublicKey &publicKey,
                                   const vector<R> &encryptedShare, const Encryption::DecryptionProof<R> &proof) {
        return this->encryptionSystem->verifyDecryption(params.encryptionParams, publicKey, encryptedShare, proof);
    }

    template<typename R>
    R PVSS<R>::reconstruct(const Params &params, const vector<vector<R> > &decryptedShares) {
        return params.prime; // TODO: Lagrange
    }
}
