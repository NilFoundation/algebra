//---------------------------------------------------------------------------//
// Copyright (c) 2020 Mikhail Komarov <nemo@nil.foundation>
// Copyright (c) 2020 Nikita Kaskov <nbering@nil.foundation>
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//---------------------------------------------------------------------------//

#ifndef ALGEBRA_FIELDS_BLS12_FR_HPP
#define ALGEBRA_FIELDS_BLS12_FR_HPP

#include <nil/algebra/fields/detail/element/fp.hpp>
#include <nil/algebra/fields/detail/params/params.hpp>

#include <nil/algebra/fields/field.hpp>

#include <nil/algebra/detail/literals.hpp>

namespace nil {
    namespace algebra {
        namespace fields {

            /*!
             * @brief
             * @tparam ModulusBits
             * @tparam GeneratorBits
             */
            template<std::size_t ModulusBits, std::size_t GeneratorBits = CHAR_BIT>
            struct bls12_scalar_field : public field<ModulusBits, GeneratorBits> { };

            template<>
            struct bls12_scalar_field<381, CHAR_BIT> : public field<381, CHAR_BIT> {
                typedef field<255, CHAR_BIT> policy_type;

                constexpr static const std::size_t modulus_bits = policy_type::modulus_bits;
                typedef typename policy_type::modulus_type modulus_type;

                constexpr static const modulus_type modulus =
                    0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001_cppui255;

                constexpr static const std::size_t generator_bits = policy_type::generator_bits;
                typedef typename policy_type::generator_type generator_type;

                constexpr static const generator_type mul_generator = 0x07;

                typedef typename detail::element_fp<detail::extension_params<bls12_scalar_field<381, CHAR_BIT>>> value_type;

                constexpr static const std::size_t arity = 1;
            };

            template<>
            struct bls12_scalar_field<377, CHAR_BIT> : public field<377, CHAR_BIT> {
                typedef field<255, CHAR_BIT> policy_type;

                constexpr static const std::size_t modulus_bits = policy_type::modulus_bits;
                typedef typename policy_type::modulus_type modulus_type;

                constexpr static const modulus_type modulus =
                    0x12AB655E9A2CA55660B44D1E5C37B00159AA76FED00000010A11800000000001_cppui253;

                constexpr static const std::size_t generator_bits = policy_type::generator_bits;
                typedef typename policy_type::generator_type generator_type;

                constexpr static const generator_type mul_generator = 0x16;

                typedef typename detail::element_fp<detail::extension_params<bls12_scalar_field<377, CHAR_BIT>>> value_type;

                constexpr static const std::size_t arity = 1;
            };

            constexpr typename bls12_scalar_field<381, CHAR_BIT>::modulus_type const bls12_scalar_field<381, CHAR_BIT>::modulus;
            constexpr typename bls12_scalar_field<377, CHAR_BIT>::modulus_type const bls12_scalar_field<377, CHAR_BIT>::modulus;

            constexpr typename bls12_scalar_field<381, CHAR_BIT>::generator_type const bls12_scalar_field<381, CHAR_BIT>::mul_generator;
            constexpr typename bls12_scalar_field<377, CHAR_BIT>::generator_type const bls12_scalar_field<377, CHAR_BIT>::mul_generator;

            template<std::size_t ModulusBits = 381, std::size_t GeneratorBits = CHAR_BIT>
            using bls12_fr = bls12_scalar_field<ModulusBits, GeneratorBits>;
            template<std::size_t ModulusBits = 377, std::size_t GeneratorBits = CHAR_BIT>
            using bls12_fr = bls12_scalar_field<ModulusBits, GeneratorBits>;

        }    // namespace fields
    }        // namespace algebra
}    // namespace nil

#endif    // ALGEBRA_FIELDS_BLS12_FR_HPP