#ifndef MY_IMPLEMENTATION_HPP
#define MY_IMPLEMENTATION_HPP
/* my_implementation.hpp
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

namespace MyEncryption {
    struct MyParams final : MyFramework::Encryption::Params {
        long l;
        long m;
        long d;
        NTL::ZZ module;
    };

    struct MyPrivateKey final : MyFramework::Encryption::PrivateKey {
        NTL::mat_ZZ trapdoor;
    };

    struct MyPublicKey final : MyFramework::Encryption::PublicKey {
        NTL::mat_ZZ A;
    };

    struct MyKeyProof final : MyFramework::Encryption::KeyProof {
        NTL::vec_ZZ x;
    };

    struct MyDecryptionProof final : MyFramework::Encryption::DecryptionProof {
        NTL::vec_ZZ e;
        NTL::mat_ZZ B;
    };

    class EncryptionType1 final : MyFramework::Encryption::EncryptionSystem {
        static void _generateTrapdoor(NTL::mat_ZZ &A, NTL::mat_ZZ &trapdoor, long n, long m, const NTL::ZZ &module,
                                      const NTL::ZZ &bound);

        static void _preSample(NTL::vec_ZZ &x, const NTL::mat_ZZ &trapdoor, const NTL::mat_ZZ &A, const NTL::vec_ZZ &b,
                               const NTL::ZZ &bound);

        static void _hash(NTL::vec_ZZ &hash, const NTL::mat_ZZ &A);

    public:
        void setup(MyFramework::Encryption::Params *&params, long securityParameter, NTL::ZZ plainBound) override;

        void generateKey(MyFramework::Encryption::KeyPair *&key, const MyFramework::Encryption::Params *params) override;

        bool verifyKey(const MyFramework::Encryption::Params *params,
                       const MyFramework::Encryption::PublicKey *publicKey,
                       const MyFramework::Encryption::KeyProof *proof) override;

        void generateEncryptionFunctionFromInput(NTL::mat_ZZ &f1, NTL::mat_ZZ &f2,
                                                 const MyFramework::Encryption::Params *params,
                                                 const MyFramework::Encryption::PublicKey *publicKey,
                                                 const NTL::vec_ZZ &plainValues,
                                                 const NTL::vec_ZZ &randomValues) override;

        void generateEncryptionFunctionFromOutput(NTL::mat_ZZ &f1, NTL::mat_ZZ &f2,
                                                  const MyFramework::Encryption::Params *params,
                                                  const MyFramework::Encryption::PublicKey *publicKey,
                                                  const NTL::vec_ZZ &cipherValues) override;

        void decrypt(MyFramework::Encryption::DecryptionProof *&proof, const MyFramework::Encryption::Params *params,
                     const MyFramework::Encryption::PublicKey *publicKey,
                     const MyFramework::Encryption::PrivateKey *privateKey, const NTL::vec_ZZ &cipherValues) override;

        bool verifyDecryption(const MyFramework::Encryption::Params *params,
                              const MyFramework::Encryption::PublicKey *publicKey, const NTL::vec_ZZ &cipherValues,
                              const MyFramework::Encryption::DecryptionProof *proof) override;
    };

    class EncryptionType2 final : MyFramework::Encryption::EncryptionSystem {
    public:
        void setup(MyFramework::Encryption::Params *&params, long securityParameter, NTL::ZZ plainBound) override;

        void generateKey(MyFramework::Encryption::KeyPair *&key, const MyFramework::Encryption::Params *params) override;

        bool verifyKey(const MyFramework::Encryption::Params *params,
                       const MyFramework::Encryption::PublicKey *publicKey,
                       const MyFramework::Encryption::KeyProof *proof) override;

        void generateEncryptionFunctionFromInput(NTL::mat_ZZ &f1, NTL::mat_ZZ &f2,
                                                 const MyFramework::Encryption::Params *params,
                                                 const MyFramework::Encryption::PublicKey *publicKey,
                                                 const NTL::vec_ZZ &plainValues,
                                                 const NTL::vec_ZZ &randomValues) override;

        void generateEncryptionFunctionFromOutput(NTL::mat_ZZ &f1, NTL::mat_ZZ &f2,
                                                  const MyFramework::Encryption::Params *params,
                                                  const MyFramework::Encryption::PublicKey *publicKey,
                                                  const NTL::vec_ZZ &cipherValues) override;

        void decrypt(MyFramework::Encryption::DecryptionProof *&proof, const MyFramework::Encryption::Params *params,
                     const MyFramework::Encryption::PublicKey *public_key,
                     const MyFramework::Encryption::PrivateKey *privateKey, const NTL::vec_ZZ &cipherValues) override;

        bool verifyDecryption(const MyFramework::Encryption::Params *params,
                              const MyFramework::Encryption::PublicKey *publicKey, const NTL::vec_ZZ &cipherValues,
                              const MyFramework::Encryption::DecryptionProof *proof) override;
    };
}

namespace MyVectorCommitment {
    struct MyParams final : MyFramework::VC::Params {
        NTL::mat_ZZ A;
        NTL::vec_ZZ v1, v2, h, h2, t;
        NTL::Mat<NTL::vec_ZZ> u;
        NTL::ZZ p, q, beta, alpha;
        long l, n;
        bool isComplement = false;
    };

    struct MyCommitment final : MyFramework::VC::Commitment {
        NTL::ZZ c;
    };

    struct MyAuxiliary final : MyFramework::VC::Auxiliary {
        NTL::vec_vec_ZZ u1, u2;
        NTL::vec_ZZ x1, x2;
    };

    struct MyOpeningProof final : MyFramework::VC::OpeningProof {
        NTL::ZZ _c;
        NTL::vec_ZZ pi, pi_eq, pi_ip;
    };

    class VectorCommitmentType1 final : MyFramework::VC::VectorCommitmentSystem {
        static void _generateTrapdoor(NTL::mat_ZZ &A, NTL::mat_ZZ &trapdoor, long n, long m, const NTL::ZZ &module,
                                      const NTL::ZZ &bound);

        static void _preSample(NTL::vec_ZZ &x, const NTL::mat_ZZ &trapdoor, const NTL::mat_ZZ &A, const NTL::vec_ZZ &b,
                               const NTL::ZZ &bound);

    public:
        void setup(MyFramework::VC::Params *&params, long securityParameter, long firstInputSize, long secondInputSize,
                   long outputSize, NTL::ZZ firstInputBound, NTL::ZZ secondInputBound,
                   NTL::ZZ coefficientBound) override;

        void commit(MyFramework::VC::Commitment *&commitment, MyFramework::VC::Auxiliary *&auxiliary,
                    const MyFramework::VC::Params *params, const NTL::vec_ZZ &firstInput,
                    const NTL::vec_ZZ &secondInput) override;

        void open(MyFramework::VC::OpeningProof *&proof, const MyFramework::VC::Params *params,
                  MyFramework::VC::Auxiliary *auxiliary, const NTL::mat_ZZ &openingFunction1,
                  const NTL::mat_ZZ &openingFunction2) override;

        bool verify(const MyFramework::VC::Params *params, const NTL::mat_ZZ &openingFunction1,
                    const NTL::mat_ZZ &openingFunction2, const MyFramework::VC::Commitment *commitment,
                    const MyFramework::VC::OpeningProof *proof) override;
    };
}

#endif //MY_IMPLEMENTATION_HPP
