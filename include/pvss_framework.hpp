#ifndef PVSS_FRAMEWORK_HPP
#define PVSS_FRAMEWORK_HPP
/* pvss_framework.hpp - a general PVSS framework
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

#include <vector>
#include <utility>
#include <memory>

namespace MyFramework {
    namespace Encryption {
        struct Params {
            int randomSize;
            int cipherTextSize;
            int randomBound;
        };

        struct PrivateKey {
        };

        struct PublicKey {
        };

        struct KeyProof {
        };

        struct KeyPair {
            PublicKey publicKey;
            PrivateKey privateKey;
            KeyProof proof;
        };

        template<typename R>
        struct DecryptionProof {
            std::vector<R> decryptedText;
        };

        template<typename R>
        class EncryptionSystem {
        public:
            virtual ~EncryptionSystem() = default;

            virtual Params setup(int securityParameter) = 0;

            virtual KeyPair generateKey(const Params &params) = 0;

            virtual bool verifyKey(const Params &params, const PublicKey &publicKey,
                                   const KeyProof &proof) = 0;

            virtual std::vector<std::vector<R> > generateEncryptionFunctionFromInput(
                const Params &params, const PublicKey &publicKey, const R &plainText,
                const std::vector<R> &randomValues) = 0;

            virtual std::vector<std::vector<R> > generateEncryptionFunctionFromOutput(
                const Params &params, const PublicKey &publicKey, const std::vector<R> &cipherText) = 0;

            virtual DecryptionProof<R> decrypt(const Params &params, const PrivateKey &privateKey,
                                               const std::vector<R> &cipherText) = 0;

            virtual bool verifyDecryption(const Params &params, const PublicKey &publicKey,
                                          const std::vector<R> &cipherText, const DecryptionProof<R> &proof) = 0;
        };
    }

    namespace VC {
        struct Params {
            int inputSize;
            int outputSize;
            int boundedInputSize;
            int bound;
        };

        struct Commitment {
        };

        struct Auxiliary {
        };

        template<typename R>
        struct OpeningProof {
            std::vector<R> output;
        };

        template<typename R>
        class VectorCommitmentSystem {
        public:
            virtual ~VectorCommitmentSystem() = default;

            virtual Params setup(int securityParameter, int inputSize, int outputSize, int boundedInputSize,
                                 int bound) = 0;

            virtual std::pair<Commitment, Auxiliary> commit(const Params &params, const std::vector<R> &input)
            = 0;

            virtual OpeningProof<R> open(const Params &params, Auxiliary &aux,
                                         const std::vector<std::vector<R> > &openingFunction) = 0;

            virtual bool verify(const Params &params, const std::vector<std::vector<R> > &openingFunction,
                                const Commitment &commitment, const OpeningProof<R> &proof) = 0;
        };
    }

    struct Params {
        int numberOfParties;
        int threshold;
        int prime;
        Encryption::Params encryptionParams;
        VC::Params vcParams;
    };

    template<typename R>
    struct DistributionProof {
        std::vector<std::vector<R> > encryptedShares;
        VC::Commitment commitment;
        VC::OpeningProof<R> proof;
    };

    template<typename R>
    class PVSS {
        std::unique_ptr<Encryption::EncryptionSystem<R> > encryptionSystem;
        std::unique_ptr<VC::VectorCommitmentSystem<R> > vectorCommitmentSystem;

    public:
        PVSS(std::unique_ptr<Encryption::EncryptionSystem<R> > encryptionSystem,
             std::unique_ptr<VC::VectorCommitmentSystem<R> > vectorCommitmentSystem) {
            this->encryptionSystem = std::move(encryptionSystem);
            this->vectorCommitmentSystem = std::move(vectorCommitmentSystem);
        }

        Params setup(int securityParameter, int numberOfParties, int threshold);

        Encryption::KeyPair generateKey(const Params &params, int index);

        bool verifyKey(const Params &params, int index, const Encryption::PublicKey &publicKey,
                       const Encryption::KeyProof &proof);

        DistributionProof<R> distribute(const Params &params,
                                        const std::vector<Encryption::PublicKey> &publicKeys, const R &secret);

        bool verifyDistribution(const Params &params, const std::vector<Encryption::PublicKey> &publicKeys,
                                const DistributionProof<R> &proof);

        Encryption::DecryptionProof<R> decryptShare(const Params &params, int index,
                                                    const Encryption::PrivateKey &privateKey,
                                                    const std::vector<R> &encryptedShare);

        bool verifyDecryption(const Params &params, int index, const Encryption::PublicKey &publicKey,
                              const std::vector<R> &encryptedShare, const Encryption::DecryptionProof<R> &proof);

        R reconstruct(const Params &params, const std::vector<std::vector<R> > &decryptedShares);
    };
}

#endif //PVSS_FRAMEWORK_HPP
