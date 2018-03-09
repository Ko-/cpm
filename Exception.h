#ifndef _EXCEPTION_
#define _EXCEPTION_

/**
 * Exception with a string expressing the reason.
 */
class Exception {
public:
    /** A string expressing the reason for the exception. */
    std::string reason;

public:
    /**
     * The default constructor.
     */
    Exception()
        : reason() {}

    /**
     * The constructor.
     * @param   aReason     A string giving the reason of the exception.
     */
    explicit Exception(const std::string& aReason)
        : reason(aReason) {}
};

#endif
