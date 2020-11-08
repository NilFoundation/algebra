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

#ifndef CRYPTO3_ALGEBRA_PAIRING_BLS12_FINAL_EXPONENTIATION_HPP
#define CRYPTO3_ALGEBRA_PAIRING_BLS12_FINAL_EXPONENTIATION_HPP

#include <nil/crypto3/algebra/pairing/detail/bls12/basic_policy.hpp>

#include <boost/multiprecision/number.hpp>
#include <boost/multiprecision/cpp_int.hpp>

namespace nil {
    namespace crypto3 {
        namespace algebra {
            namespace pairing {
                namespace policies {
                    namespace detail {

                        /*************************  FINAL EXPONENTIATIONS  ***********************************/

                        template<std::size_t ModulusBits = 381>
                        class bls12_final_exponentiation_basic_functions;

                        template<>
                        class bls12_final_exponentiation_basic_functions<381>{
                            constexpr static const std::size_t modulus_bits = 381;
                            using policy_type = pairing::detail::bls12_basic_policy<modulus_bits>;

                        public:
                            using gt = typename policy_type::gt;

                            static gt final_exponentiation_first_chunk(const gt &elt) {

                                /*
                                  Computes result = elt^((q^6-1)*(q^2+1)).
                                  Follows, e.g., Beuchat et al page 9, by computing result as follows:
                                     elt^((q^6-1)*(q^2+1)) = (conj(elt) * elt^(-1))^(q^2+1)
                                  More precisely:
                                  A = conj(elt)
                                  B = elt.inversed()
                                  C = A * B
                                  D = C.Frobenius_map(2)
                                  result = D * C
                                */

                                const gt A = elt.unitary_inversed();
                                const gt B = elt.inversed();
                                const gt C = A * B;
                                const gt D = C.Frobenius_map(2);
                                const gt result = D * C;

                                return result;
                            }

                            static gt exp_by_z(const gt &elt) {

                                gt result = elt.cyclotomic_exp(policy_type::final_exponent_z);
                                if (policy_type::final_exponent_is_z_neg) {
                                    result = result.unitary_inversed();
                                }

                                return result;
                            }

                            static gt final_exponentiation_last_chunk(const gt &elt) {

                                const gt A = elt.cyclotomic_squared();    // elt^2
                                const gt B = A.unitary_inversed();        // elt^(-2)
                                const gt C = exp_by_z(elt);               // elt^z
                                const gt D = C.cyclotomic_squared();      // elt^(2z)
                                const gt E = B * C;                       // elt^(z-2)
                                const gt F = exp_by_z(E);                 // elt^(z^2-2z)
                                const gt G = exp_by_z(F);                 // elt^(z^3-2z^2)
                                const gt H = exp_by_z(G);                 // elt^(z^4-2z^3)
                                const gt I = H * D;                       // elt^(z^4-2z^3+2z)
                                const gt J = exp_by_z(I);                 // elt^(z^5-2z^4+2z^2)
                                const gt K = E.unitary_inversed();        // elt^(-z+2)
                                const gt L = K * J;                       // elt^(z^5-2z^4+2z^2) * elt^(-z+2)
                                const gt M = elt * L;                     // elt^(z^5-2z^4+2z^2) * elt^(-z+2) * elt
                                const gt N = elt.unitary_inversed();      // elt^(-1)
                                const gt O = F * elt;                     // elt^(z^2-2z) * elt
                                const gt P = O.Frobenius_map(3);          // (elt^(z^2-2z) * elt)^(q^3)
                                const gt Q = I * N;                       // elt^(z^4-2z^3+2z) * elt^(-1)
                                const gt R = Q.Frobenius_map(1);          // (elt^(z^4-2z^3+2z) * elt^(-1))^q
                                const gt S = C * G;                       // elt^(z^3-2z^2) * elt^z
                                const gt T = S.Frobenius_map(2);          // (elt^(z^3-2z^2) * elt^z)^(q^2)
                                const gt U = T * P;    // (elt^(z^2-2z) * elt)^(q^3) * (elt^(z^3-2z^2) * elt^z)^(q^2)
                                const gt V = U * R;    // (elt^(z^2-2z) * elt)^(q^3) * (elt^(z^3-2z^2) * elt^z)^(q^2) *
                                                       // (elt^(z^4-2z^3+2z) * elt^(-1))^q
                                const gt W =
                                    V * M;    // (elt^(z^2-2z) * elt)^(q^3) * (elt^(z^3-2z^2) * elt^z)^(q^2) *
                                              // (elt^(z^4-2z^3+2z) * elt^(-1))^q * elt^(z^5-2z^4+2z^2) * elt^(-z+2) * elt

                                return W;
                            }
                        };
                    }    // namespace detail


                    template<std::size_t ModulusBits = 381>
                    class bls12_final_exponentiation;

                    template<>
                    class bls12_final_exponentiation<381>{
                        constexpr static const std::size_t modulus_bits = 381;
                        using policy_type = pairing::detail::bls12_basic_policy<modulus_bits>;

                    public:
                        using gt = typename policy_type::gt;
                        
                        gt operator()(const gt &elt) {
                            /* OLD naive version:
                                gt result = elt^final_exponent;
                            */
                            gt A = detail::bls12_final_exponentiation_basic_functions<modulus_bits>::final_exponentiation_first_chunk(elt);
                            gt result = detail::bls12_final_exponentiation_basic_functions<modulus_bits>::final_exponentiation_last_chunk(A);

                            return result;
                        }
                    };
                }    // namespace policies
            }        // namespace pairing
        }            // namespace algebra
    }                // namespace crypto3
}    // namespace nil
#endif    // CRYPTO3_ALGEBRA_PAIRING_BLS12_FINAL_EXPONENTIATION_HPP
