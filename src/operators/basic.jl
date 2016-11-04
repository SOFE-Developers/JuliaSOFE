# EXPORTS
export ScalarOperator, VectorOperator, MatrixOperator

abstract ScalarOperator <: AbstractOperator
abstract VectorOperator <: AbstractOperator
abstract MatrixOperator <: AbstractOperator

#-------------------------#
# Identity Operator Types #
#-------------------------#
type id <: ScalarOperator end
type Id <: VectorOperator end
type ID <: MatrixOperator end

#-------------------------#
# Gradient Operator Types #
#-------------------------#
type grad <: ScalarOperator end
type Grad <: VectorOperator end
type GRAD <: MatrixOperator end

#---------------------------#
# Divergence Operator Types #
#---------------------------#
type div <: ScalarOperator end
type Div <: VectorOperator end
type DIV <: MatrixOperator end
