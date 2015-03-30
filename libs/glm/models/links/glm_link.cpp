#include <glm/models/links/glm_link.hpp>

#include <glm/models/links/identity.hpp>
#include <glm/models/links/log.hpp>
#include <glm/models/links/logc.hpp>
#include <glm/models/links/odds.hpp>
#include <glm/models/links/logit.hpp>

glm_link *
make_link(const std::string &link_name)
{
    if( link_name == "identity" )
    {
        return new identity_link( );
    }
    else if( link_name == "log" )
    {
        return new log_link( );
    }
    else if( link_name == "logc" )
    {
        return new logc_link( );
    }
    else if( link_name == "odds" )
    {
        return new odds_link( );
    }
    else if( link_name == "logit" )
    {
        return new logit_link( );
    }
    else
    {
        return NULL;
    }
}
