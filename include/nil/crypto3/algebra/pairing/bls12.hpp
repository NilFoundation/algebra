//---------------------------------------------------------------------------//
// Copyright (c) 2020 Mikhail Komarov <nemo@nil.foundation>
// Copyright (c) 2020 Nikita Kaskov <nbering@nil.foundation>
//
// MIT License
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//---------------------------------------------------------------------------//

#ifndef CRYPTO3_ALGEBRA_PAIRING_BLS12_POLICY_HPP
#define CRYPTO3_ALGEBRA_PAIRING_BLS12_POLICY_HPP

#include <nil/crypto3/algebra/pairing/detail/bls12/basic_policy.hpp>
#include <nil/crypto3/algebra/pairing/policies/bls12/final_exponentiation.hpp>

namespace nil {
    namespace crypto3 {
        namespace algebra {
            namespace curves {

                template<std::size_t ModulusBits>
                struct bls12;

            }    // namespace curves
            namespace pairing {

                

                template<std::size_t ModulusBits, typename PairingFunctions, 
                    typename FinalExponentiation = policies::bls12_final_exponentiation<ModulusBits>>
                struct bls12_pairing_policy {
                    using policy_type = detail::bls12_basic_policy<ModulusBits>;
                    using functions_policy = PairingFunctions;
                public:
                    using number_type = typename policy_type::number_type;

                    constexpr static const typename policy_type::number_type pairing_loop_count =
                        policy_type::ate_loop_count;

                    using Fp_type = typename policy_type::Fp_field;
                    using G1_type = typename policy_type::g1;
                    using G2_type = typename policy_type::g2;
                    using Fq_type = typename policy_type::Fq_field;
                    using Fqe_type = typename policy_type::Fqe_field;
                    using Fqk_type = typename policy_type::Fqk_field;
                    using GT_type = typename policy_type::gt;

                    using G1_precomp = typename functions_policy::g1_precomp;
                    using G2_precomp = typename functions_policy::g2_precomp;

                    static inline typename functions_policy::g1_precomp precompute_g1(const typename policy_type::g1 &P) {
                        return functions_policy::precompute_g1(P);
                    }

                    static inline typename functions_policy::g2_precomp precompute_g2(const typename policy_type::g2 &Q) {
                        return functions_policy::precompute_g2(Q);
                    }

                    static inline typename policy_type::gt pairing(const typename policy_type::g1 &P,
                                                                   const typename policy_type::g2 &Q) {
                        return functions_policy::pairing(P, Q);
                    }

                    static inline typename policy_type::gt reduced_pairing(const typename policy_type::g1 &P,
                                                                           const typename policy_type::g2 &Q) {
                        return functions_policy::reduced_pairing(P, Q);
                    }

                    static inline typename policy_type::gt
                        double_miller_loop(const typename functions_policy::g1_precomp &prec_P1,
                                           const typename functions_policy::g2_precomp &prec_Q1,
                                           const typename functions_policy::g1_precomp &prec_P2,
                                           const typename functions_policy::g2_precomp &prec_Q2) {
                        return functions_policy::double_miller_loop(prec_P1, prec_Q1, prec_P2, prec_Q2);
                    }

                    using final_exponentiation = FinalExponentiation;

                    static inline typename policy_type::gt miller_loop(const typename functions_policy::g1_precomp &prec_P,
                                                                       const typename functions_policy::g2_precomp &prec_Q) {
                        return functions_policy::miller_loop(prec_P, prec_Q);
                    }
                };

                template<std::size_t ModulusBits, typename PairingFunctions, typename FinalExponentiation>
                constexpr typename bls12_pairing_policy<ModulusBits, PairingFunctions, FinalExponentiation>::number_type const
                    bls12_pairing_policy<ModulusBits, PairingFunctions, FinalExponentiation>::pairing_loop_count;
            }    // namespace pairing
        }        // namespace algebra
    }            // namespace crypto3
}    // namespace nil
#endif    // CRYPTO3_ALGEBRA_PAIRING_BLS12_POLICY_HPP