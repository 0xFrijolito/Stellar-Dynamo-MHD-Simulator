#ifndef VECTOR_H
#define VECTOR_H

#include "core/defines.h"

/**
 * @file vector.h
 * @brief Operaciones con vectores tridimensionales abstractos (Vector3).
 *
 * Este módulo define el tipo Vector3 y proporciona funciones para operar con
 * vectores en un espacio tridimensional general. No se asume ninguna métrica
 * específica ni un sistema de coordenadas fijo; el significado de las
 * componentes queda determinado por el contexto donde se utiliza el vector.
 */

    /**
     * @struct Vector3
     * @brief Representa un vector tridimensional abstracto.
     *
     * Cada componente del vector es un número de punto flotante de doble precisión.
     * No se asume una correspondencia con coordenadas cartesianas (x, y, z), ya que
     * el significado de cada componente depende de la geometría y métrica del espacio
     * definido por quien utilice esta estructura.
     */
    typedef struct Vector3 {
        f64 x1; /**< Primera componente del vector. */
        f64 x2; /**< Segunda componente del vector. */
        f64 x3; /**< Tercera componente del vector. */
    } Vector3;

/**
 * @brief Genera un vector tridimensional con componentes aleatorias.
 *
 * La distribución de los valores aleatorios puede depender de la implementación.
 *
 * @return Un Vector3 con valores aleatorios en sus tres componentes.
 */
Vector3 vector3_random();

/**
 * @brief Crea un vector tridimensional con todas las componentes en cero.
 *
 * @return Un Vector3 cuyas componentes son todas 0.
 */
Vector3 vector3_zeros();

/**
 * @brief Crea un vector tridimensional con todas las componentes iguales a un valor dado.
 *
 * @param value Valor constante a asignar a cada componente del vector.
 * @return Un Vector3 con todas sus componentes iguales a `value`.
 */
Vector3 vector3_fill(f64 value);

/**
 * @brief Suma componente a componente dos vectores tridimensionales.
 *
 * @param v1 Primer vector.
 * @param v2 Segundo vector.
 * @return Vector resultante de sumar `v1` y `v2`.
 */
Vector3 vector3_add(Vector3 v1, Vector3 v2);

/**
 * @brief Resta componente a componente dos vectores tridimensionales.
 *
 * @param v1 Primer vector.
 * @param v2 Segundo vector.
 * @return Vector resultante de `v1 - v2`.
 */
Vector3 vector3_sub(Vector3 v1, Vector3 v2);

/**
 * @brief Multiplica cada componente de un vector por un escalar.
 *
 * @param vector Vector a escalar.
 * @param value Escalar por el que se multiplican las componentes.
 * @return Vector resultante de la multiplicación.
 */
Vector3 vector3_scalar_mul(Vector3 vector, f64 value);

/**
 * @brief Calcula el producto punto entre dos vectores.
 *
 * El resultado se almacena en una variable externa pasada por referencia.
 * La métrica que define este producto depende del contexto donde se use.
 *
 * @param v1 Primer vector.
 * @param v2 Segundo vector.
 * @param output Variable donde se almacenará el resultado del producto punto.
 */
void vector3_dot_product(Vector3 v1, Vector3 v2, f32 output);

/**
 * @brief Calcula el producto cruz entre dos vectores.
 *
 * El resultado se almacena en un vector externo pasado por referencia. Su definición
 * puede depender del álgebra del espacio en el que se trabaja.
 *
 * @param v1 Primer vector.
 * @param v2 Segundo vector.
 * @param output Vector donde se almacenará el resultado del producto cruz.
 */
void vector3_cross_product(Vector3 v1, Vector3 v2, Vector3 output);

/**
 * @brief Copia profundamente un vector en otro.
 *
 * Copia el contenido del vector fuente `v2` al vector destino `v1`.
 *
 * @param v1 Puntero al vector destino.
 * @param v2 Puntero al vector fuente.
 */
void vector3_deep_copy(Vector3 *v1, Vector3 *v2);

#endif
