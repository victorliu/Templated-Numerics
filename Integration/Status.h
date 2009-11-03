#ifndef _INTEGRATION_STATUS_H_
#define _INTEGRATION_STATUS_H_

// Preprocessor macros

namespace Integration{

template <typename NumericType>
class StatusCallback{
public:
	struct Status{
		size_t num_evals;
		NumericType current_value;
		NumericType current_error;
	};
	virtual void operator()(const Status &status) = 0;
	virtual size_t NumUpdates() const{ return 100; }
	virtual bool WantUpdates() const{ return true; }
};

template <typename NumericType>
class DefaultStatusCallback : StatusCallback<NumericType>{
public:
	void operator()(const Status &status){}
	bool WantUpdates() const{ return false; }
};

}; // namespace Integration

#endif // _INTEGRATION_STATUS_H_
