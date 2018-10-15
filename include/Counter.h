#ifndef COUNTER_H
#define COUNTER_H
template <typename T>
struct Counter
{
	Counter()
	{
		objects_created++;
		objects_alive++;
	}

	virtual ~Counter()
	{
		--objects_alive;
	}
	static int objects_created;
	static int objects_alive;
};
template <typename T> int Counter<T>::objects_created(0);
template <typename T> int Counter<T>::objects_alive(0);

#endif