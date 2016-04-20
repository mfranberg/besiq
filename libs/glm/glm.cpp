#include <glm/glm.hpp>

#include <glm/lm.hpp>
#include <glm/irls.hpp>

arma::vec
glm_fit(const arma::mat &X, const arma::vec &y, const arma::uvec &missing, const glm_model &model, glm_info &output, bool fast_inversion)
{
    if( model.get_name( ) == "normal" && model.get_link( ).get_name( ) == "identity" )
    {
        return lm( X, y, missing, model, output );
    }
    else
    {
        return irls( X, y, missing, model, output, fast_inversion );
    }
}
