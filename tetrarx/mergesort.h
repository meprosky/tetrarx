#pragma once

//взято с ru.wikibooks.org/wiki/
// Слияние двух упорядоченных массивов
// m1 - Первый массив
// m2 - Второй массив
// l1 - Длина первого массива
// l2 - Длина второго массива
// Возвращает объединённый массив
              
template <class T> T* merge(T *m1, T* m2, int l1, int l2)
{
    T* ret = new T[l1+l2];
    int n = 0;
    // Сливаем массивы, пока один не закончится
    while (l1 && l2)
    {
        if (*m1 < *m2)
        {
            ret[n] = *m1;
            m1++;
            l1--;
        }
        else
        {
            ret[n] = *m2;
            m2++;
            l2--;
        }
        n++;
    }
    // Если закончился первый массив
    if (l1 == 0)
    {
        for (int i=0; i<l2; i++)
        {
            ret[n++] = *m2++;
        }
    }
    // Если закончился второй массив
    else
    {
        for (int i=0; i<l1; i++)
        {
            ret[n++] = *m1++;
        }
    }
    return ret;
}
 
// Функция восходящего слияния
template <class T> void mergeSort(T * mas, int len)
{
    int n = 1, l, ost;
    T * mas1;
    while (n < len)
    {
        l = 0;
        while (l < len)
        {
            if (l+n >= len) break;
            ost = (l+n*2>len) ? (len-(l+n)) : n;
            mas1 = merge(mas + l, mas + l + n, n, ost);
            for (int i = 0; i<n+ost; i++) mas[l+i] = mas1[i];
            delete [] mas1;
            l += n*2;
        }
        n *= 2;
    }
}


