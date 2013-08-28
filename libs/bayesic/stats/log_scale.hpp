#ifndef __LOG_SCALE_H__
#define __LOG_SCALE_H__

#include <cmath>
#include <cfloat>
#include <stdexcept>

/**
 * This class represents a value on the log scale that supports
 * the common arithmetic operators.
 */
template <class T>
class log_scale
{
public:
    /**
     * Constructor.
     *
     * @param value The value to log.
     */
    log_scale(T value)
    {
        if( value != 0.0 )
        {
            m_log_value = log( value );
            m_is_zero = false;
        }
        else
        {
            m_log_value = -DBL_MAX;
            m_is_zero = true;
        }
    }

    /**
     * Secondary constructor.
     *
     * Intializes the class directly from a value already on the
     * log scale.
     *
     * @param log_value A already logged value.
     *
     * @return The initialized object.
     */
    template <class K>
    static log_scale<K> from_log(K log_value)
    {
        log_scale<K> from_log( 1.0 );
        from_log.m_log_value = log_value;

        return from_log;
    }

    /**
     * Assignment operator.
     */
    const log_scale<T> &operator=(const log_scale<T> &other)
    {
        if( this != &other )
        {
            m_log_value = other.m_log_value;
            m_is_zero = other.m_is_zero;
        }
        return *this;
    }

    /* Arithmetic operators +, -, * and /. */
    const log_scale<T> &operator+=(const log_scale<T> &other)
    {
        if( !m_is_zero && !other.m_is_zero )
        {
            T diff = other.m_log_value - m_log_value;
            m_log_value += log( 1.0 + exp( diff ) );
        }
        else if( m_is_zero )
        {
            m_log_value = other.m_log_value;
            m_is_zero = other.m_is_zero;
        }

        return *this;
    }

    const log_scale<T> &operator-=(const log_scale<T> &other)
    {
        if( !m_is_zero && !other.m_is_zero )
        {
            T diff = m_log_value - other.m_log_value;
            if( diff < 0.0 )
            {
                throw std::domain_error( "log_scale: Cannot represent negative values on log scale." );
            }
            m_log_value = other.m_log_value + log( exp( diff ) - 1.0 );
        }
        else if( m_is_zero && !other.m_is_zero )
        {
            throw std::domain_error( "log_scale: Cannot represent negative values on log scale." );
        }

        return *this;
    }

    const log_scale<T> &operator*=(const log_scale<T> &other)
    {
        if( !m_is_zero && !other.m_is_zero )
        {
            m_log_value += other.m_log_value;
        }
        else
        {
            m_log_value = -DBL_MAX;
            m_is_zero = true;
        }

        return *this;
    }

    const log_scale<T> &operator/=(const log_scale<T> &other)
    {
        /* Note: Order matters */
        if( !m_is_zero && !other.m_is_zero )
        {
            m_log_value -= other.m_log_value;
        }
        else if( other.m_is_zero )
        {
            throw std::domain_error( "log_scale: Division with zero." );
        }
        else
        {
            m_log_value = -DBL_MAX;
            m_is_zero = true;
        }

        return *this;
    }

    /**
     * Returns the un-logged value.
     * 
     * @return the un-logged value.
     */
    T value() const
    {
        if( !m_is_zero )
        {
            return exp( m_log_value );
        }
        else
        {
            return 0.0;
        }
    }

    /**
     * Returns the logged value.
     * 
     * @return the logged value.
     */
    T log_value() const
    {
        if( !m_is_zero )
        {
            return m_log_value;
        }
        else
        {
            return -DBL_MAX;
        }
    }

private:
    /**
     * The underlying logged value.
     */
    T m_log_value;

    /**
     * Determines whether this should be zero.
     */
    bool m_is_zero;
};

/* Binaray arithmetic operators for +, -, * and /. */
template <class T>
log_scale<T> operator+(const log_scale<T> &lhs, const log_scale<T> &rhs)
{
    return log_scale<T>( lhs ) += rhs;
}

template <class T>
log_scale<T> operator+(const T &lhs, const log_scale<T> &rhs)
{
    return log_scale<T>( lhs ) += rhs;
}

template <class T>
log_scale<T> operator+(const log_scale<T> &lhs, const T &rhs)
{
    return log_scale<T>( lhs ) += log_scale<T>( rhs );
}

template <class T>
log_scale<T> operator-(const log_scale<T> &lhs, const log_scale<T> &rhs)
{
    return log_scale<T>( lhs ) -= rhs;
}

template <class T>
log_scale<T> operator-(const T &lhs, const log_scale<T> &rhs)
{
    return log_scale<T>( lhs ) -= rhs;
}

template <class T>
log_scale<T> operator-(const log_scale<T> &lhs, const T &rhs)
{
    return log_scale<T>( lhs ) -= log_scale<T>( rhs );
}


template <class T>
log_scale<T> operator/(const log_scale<T> &lhs, const log_scale<T> &rhs)
{
    return log_scale<T>( lhs ) /= rhs;
}

template <class T>
log_scale<T> operator/(const T &lhs, const log_scale<T> &rhs)
{
    return log_scale<T>( lhs ) /= rhs;
}

template <class T>
log_scale<T> operator/(const log_scale<T> &lhs, const T &rhs)
{
    return log_scale<T>( lhs ) /= log_scale<T>( rhs );
}

template <class T>
log_scale<T> operator*(const log_scale<T> &lhs, const log_scale<T> &rhs)
{
    return log_scale<T>( lhs ) *= rhs;
}

template <class T>
log_scale<T> operator*(const T &lhs, const log_scale<T> &rhs)
{
    return log_scale<T>( lhs ) *= rhs;
}

template <class T>
log_scale<T> operator*(const log_scale<T> &lhs, const T &rhs)
{
    return log_scale<T>( lhs ) *= log_scale<T>( rhs );
}

/* Comparisson operators, not sure about how do this nicely
   without boost::less_than_comparable. */
template <class T>
bool operator<(const log_scale<T> &lhs, const log_scale<T> &rhs)
{
    return lhs.log_value( ) < rhs.log_value( );
}

template <class T>
bool operator<(const T &lhs, const log_scale<T> &rhs)
{
    return log_scale<T>( lhs ) < rhs;
}

template <class T>
bool operator<(const log_scale<T> &lhs, const T &rhs)
{
    return lhs < log_scale<T>( rhs );
}

template <class T>
bool operator<=(const log_scale<T> &lhs, const log_scale<T> &rhs)
{
    return lhs.log_value( ) <= rhs.log_value( );
}

template <class T>
bool operator<=(const T &lhs, const log_scale<T> &rhs)
{
    return log_scale<T>( lhs ) <= rhs;
}

template <class T>
bool operator<=(const log_scale<T> &lhs, const T &rhs)
{
    return lhs <= log_scale<T>( rhs );
}

template <class T>
bool operator>(const log_scale<T> &lhs, const log_scale<T> &rhs)
{
    return lhs.log_value( ) > rhs.log_value( );
}

template <class T>
bool operator>(const T &lhs, const log_scale<T> &rhs)
{
    return log_scale<T>( lhs ) > rhs;
}

template <class T>
bool operator>(const log_scale<T> &lhs, const T &rhs)
{
    return lhs > log_scale<T>( rhs );
}

template <class T>
bool operator>=(const log_scale<T> &lhs, const log_scale<T> &rhs)
{
    return lhs.log_value( ) >= rhs.log_value( );
}

template <class T>
bool operator>=(const T &lhs, const log_scale<T> &rhs)
{
    return log_scale<T>( lhs ) >= rhs;
}

template <class T>
bool operator>=(const log_scale<T> &lhs, const T &rhs)
{
    return lhs >= log_scale<T>( rhs );
}

template <class T>
bool operator==(const log_scale<T> &lhs, const log_scale<T> &rhs)
{
    return lhs.log_value( ) == rhs.log_value( );
}

template <class T>
bool operator==(const T &lhs, const log_scale<T> &rhs)
{
    return log_scale<T>( lhs ) == rhs;
}

template <class T>
bool operator==(const log_scale<T> &lhs, const T &rhs)
{
    return lhs == log_scale<T>( rhs );
}

template <class T>
bool operator!=(const log_scale<T> &lhs, const log_scale<T> &rhs)
{
    return lhs.log_value( ) != rhs.log_value( );
}

template <class T>
bool operator!=(const T &lhs, const log_scale<T> &rhs)
{
    return log_scale<T>( lhs ) != rhs;
}

template <class T>
bool operator!=(const log_scale<T> &lhs, const T &rhs)
{
    return lhs != log_scale<T>( rhs );
}

/* Convenience definitions for common types. */
typedef log_scale<double> log_double;
typedef log_scale<float> log_float;
typedef log_scale<int> log_int;

#endif /* End of __LOG_SCALE_H__ */
