// ReSharper disable CppDoxygenUndocumentedParameter
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
#include <memory>
#include <NTL/mat_ZZ.h>

namespace MyFramework {
    namespace Encryption {
        struct Params {
            long plainSize;
            NTL::ZZ plainBound;
            long randomSize;
            NTL::ZZ randomBound;
            NTL::ZZ coefficientBound;
            long cipherSize;
        };

        struct PrivateKey {
        };

        struct PublicKey {
        };

        struct KeyProof {
        };

        struct KeyPair {
            PublicKey *publicKey;
            PrivateKey *privateKey;
            KeyProof *proof;
        };

        struct DecryptionProof {
            NTL::vec_ZZ decryptedValues;
        };

        struct EncryptionSystem {
            virtual ~EncryptionSystem() = default;

            virtual void setup(const Params *params, const long securityParameter, const NTL::ZZ plainBound) = 0;

            virtual void generateKey(const KeyPair *key, const Params *params) = 0;

            virtual bool verifyKey(const Params *params, const PublicKey *publicKey, const KeyProof *proof) = 0;

            /**
             * @return @param f1 Function coefficients for plain value
             * @return @param f2 Function coefficients for random value
             *
             * @details encryptedValues = f1 * plainValues + f2 * randomValues
             */
            virtual void generateEncryptionFunctionFromInput(NTL::mat_ZZ &f1, NTL::mat_ZZ &f2,
                                                             const Params *params,
                                                             const PublicKey *publicKey,
                                                             const NTL::vec_ZZ &plainValues,
                                                             const NTL::vec_ZZ &randomValues) = 0;

            virtual void generateEncryptionFunctionFromOutput(NTL::mat_ZZ &f1, NTL::mat_ZZ &f2,
                                                              const Params *params,
                                                              const PublicKey *publicKey,
                                                              const NTL::vec_ZZ &cipherValues) = 0;

            virtual void decrypt(const DecryptionProof *proof, const Params *params, const PublicKey *publicKey,
                                 const PrivateKey *privateKey, const NTL::vec_ZZ &cipherValues) = 0;

            virtual bool verifyDecryption(const Params *params, const PublicKey *publicKey,
                                          const NTL::vec_ZZ &cipherValues, const DecryptionProof *proof) = 0;
        };
    }

    namespace VC {
        struct Params {
            long firstInputSize;
            long secondInputSize;
            long secondInputPartitionSize;
            long outputSize;
            NTL::ZZ firstInputBound;
            NTL::ZZ secondInputBound;
            NTL::ZZ coefficientBound;
        };

        struct Commitment {
        };

        struct Auxiliary {
        };

        struct OpeningProof {
            NTL::vec_ZZ output;
        };

        struct VectorCommitmentSystem {
            virtual ~VectorCommitmentSystem() = default;

            virtual void setup(const Params *params, const long securityParameter, const long firstInputSize,
                               const long secondInputSize, const long outputSize, const NTL::ZZ firstInputBound,
                               const NTL::ZZ secondInputBound, const NTL::ZZ coefficientBound) = 0;

            virtual void commit(const Commitment *commitment, const Auxiliary *auxiliary, const Params *params,
                                const NTL::vec_ZZ &firstInput, const NTL::vec_ZZ &secondInput) = 0;

            /**
             * @param openingFunction1 Function coefficients for first input
             * @param openingFunction2 Function coefficients for second input
             *
             * @details output = openingFunction1 * firstInput + openingFunction2 * secondInput
             */
            virtual void open(const OpeningProof *proof, const Params *params, const Auxiliary *auxiliary,
                              const NTL::mat_ZZ &openingFunction1, const NTL::mat_ZZ &openingFunction2) = 0;

            virtual bool verify(const Params *params, const NTL::mat_ZZ &openingFunction1,
                                const NTL::mat_ZZ &openingFunction2, const Commitment *commitment,
                                const OpeningProof *proof) = 0;
        };
    }

    struct Params {
        long numberOfParties;
        long threshold;
        NTL::ZZ prime;
        Encryption::Params *encryptionParams;
        VC::Params *vcParams;
    };

    struct DistributionProof {
        std::vector<NTL::vec_ZZ> encryptedShares;
        VC::Commitment *commitment;
        VC::OpeningProof *proof;
    };

    struct DecryptionProof {
        Encryption::DecryptionProof *proof;
        NTL::ZZ decryptedShare;
    };

    template<typename EncryptionType, typename VectorCommitmentType>
    class PVSS final {
        static_assert(std::is_base_of_v<Encryption::EncryptionSystem, EncryptionType>,
                      "EncryptionType must inherit MyFramework::Encryption::EncryptionSystem");
        static_assert(std::is_base_of_v<VC::VectorCommitmentSystem, VectorCommitmentType>,
                      "VectorCommitmentType must inherit MyFramework::VC::VectorCommitmentSystem");

        std::unique_ptr<EncryptionType> encryptionSystem;
        std::unique_ptr<VectorCommitmentType> vectorCommitmentSystem;

    public:
        PVSS(std::unique_ptr<EncryptionType> encryptionSystem,
             std::unique_ptr<VectorCommitmentType> vectorCommitmentSystem): encryptionSystem(
                                                                                std::move(encryptionSystem)),
                                                                            vectorCommitmentSystem(
                                                                                std::move(vectorCommitmentSystem)) {
        }

        void setup(Params &params, const long securityParameter, const long numberOfParties,
                   const long threshold);

        void generateKey(Encryption::KeyPair &key, const Params &params);

        bool verifyKey(const Params &params, const Encryption::PublicKey *publicKey, const Encryption::KeyProof *proof);

        void distribute(DistributionProof &proof, const Params &params,
                        const std::vector<Encryption::PublicKey *> &publicKeys, const NTL::ZZ &secret);

        bool verifyDistribution(const Params &params, const std::vector<Encryption::PublicKey *> &publicKeys,
                                const DistributionProof &proof);

        void decryptShare(DecryptionProof &proof, const Params &params, const Encryption::PublicKey *publicKey,
                          const Encryption::PrivateKey *privateKey, const NTL::vec_ZZ &encryptedShare);

        bool verifyDecryption(const Params &params, const Encryption::PublicKey *publicKey,
                              const NTL::vec_ZZ &encryptedShare, const DecryptionProof &proof);

        void reconstruct(NTL::ZZ &reconstruction, const Params &params, const std::vector<NTL::ZZ> &decryptedShares);
    };
}

#endif //PVSS_FRAMEWORK_HPP
